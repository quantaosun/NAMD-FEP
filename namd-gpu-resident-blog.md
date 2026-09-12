# NAMD GPU-resident: a 41× FEP speedup, and why more CPU threads made it *slower*

**Or: the day I learned throughput can be inversely proportional to core count.**

*TL;DR — a NAMD 3.0.3 free-energy run that was taking ~25–30 days now finishes in
~14 hours. Half the fix was a single compile flag I'd been missing. The other half
was realizing that, in GPU-resident mode, adding CPU threads made the simulation
run *slower*, almost exactly proportionally.*

---

## The mystery

I had a 64,651-atom protein–ligand complex on a Tesla V100, running NAMD 3.0.3
alchemical free-energy (FEP) calculations. It was doing **~1.7 ns/day** plain MD
and **~1.2 ns/day** FEP. The GPU sat at **~2% utilization and 534 MB**. The CPU
sat at ~9%. Nothing was busy, and nothing was fast.

At that rate a 16-window × 2 ns production run would take **25–30 days**. That's
not a benchmark — that's a career.

Whenever I tried the obvious fix, NAMD told me:

```
FATAL ERROR: GPUresident not supported on regular multicore builds
```

---

## Fix #1: the flag that compiles the feature out

That error is the whole story in disguise. `GPUresident` (and its deprecated
spelling `CUDASOAintegrate`) is gated at **compile time** by a single macro,
`NODEGROUP_FORCE_REGISTER`:

```c
// src/SimParameters.C:5008
#ifndef NODEGROUP_FORCE_REGISTER
  NAMD_die("GPUresident not supported on regular multicore builds");
#endif
```

And that macro is emitted by **exactly one** configure option:

```bash
# config:1052 — only under `use_sn` (single-node)
if (( $use_sn )); then
    echo "EXTRADEFINES += -DNODEGROUP_FORCE_REGISTER "
    echo "CUDAFLAGS     += -DNODEGROUP_FORCE_REGISTER"
fi
```

I had built NAMD with `--with-cuda`. What I needed was `--with-single-node-cuda`.
The difference is invisible unless you go looking — but `--with-cuda` compiles the
GPU-resident path **out of the binary entirely**. The GPU was never going to
accelerate anything; it could only sit there offloading the occasional kernel.

The fix was a rebuild, not a rewrite:

```bash
./config Linux-x86_64-g++.gpuresident \
  --charm-arch multicore-linux-x86_64 \
  --with-single-node-cuda \
  --cuda-prefix /usr/local/cuda-11.8 \
  --cuda-gencode arch=compute_70,code=sm_70
make -j32   # ~13 minutes on 64 cores
```

That alone took the same FEP window from **2.26 → 8.73 ns/day** (~3.9×).

> A note I'd written to myself earlier blamed the Charm++ architecture and proposed
> rebuilding with `verbs`/`mpi`. That was wrong. `multicore-linux-x86_64` is the
> *correct* arch for single-node GPU-resident; the existing Charm++ build was reused
> untouched. Check your assumptions before you rebuild Charm++.

---

## Fix #2: the counterintuitive one — *fewer* threads is faster

This is the part I want people to steal.

Out of habit (and because it's correct for the *old* offload build), I ran with
`+p8`. In GPU-resident mode, essentially **all** of the force computation lives on
the GPU. The CPU's job is mostly to orchestrate — and extra CPU "processes" don't
parallelize that orchestration, they just **contend on synchronization**.

Measured on the same 64,651-atom complex, same config, only the PE count changed:

| threads | ns/day |
|---------|--------|
| **`+p1`** | **75.79** |
| `+p2` | 59.54 |
| `+p4` | 27.22 |
| `+p8` | 13.22 |
| `+p16` | 6.29 |
| `+p32` | 2.70 |

It's almost exactly inversely proportional. `+p2` is ~½ of `+p1`. `+p4` is ~½ of
`+p2`. `+p8` is ~½ of `+p4`. Every time you double the thread count, you roughly
halve the throughput. The `+p8` habit was costing **~5.7×**.

For fun I tried the "tuning" knobs people recommend for GPU-resident —
`stepspercycle 400`, `pairlistsPerCycle 40`, `margin 4.0`. Result: **noise**
(75.66 vs 75.79 ns/day). PE count was the entire story.

**The rule I took away: in GPU-resident mode, one PE, one GPU. `+p1`.**

---

## The result

| | old build, `+p8` | GPU-resident, `+p1` | speedup |
|---|---|---|---|
| plain MD | 3.57 ns/day | **75.79 ns/day** | **21×** |
| FEP (single window) | 2.26 ns/day | **49.65 ns/day** | **22×** |

Versus the *original* broken baseline (~1.2 ns/day FEP), that's **~41×**.

Production numbers, measured on the real generated configs rather than the
benchmark:

- **complex** (64,651 atoms): **47 ns/day**
- **solvent** (5,013 atoms): **59 ns/day** — only 13× smaller but **not** 13× faster;
  small systems on the GPU are latency-bound, so it saturates around ~59 ns/day.

The full job — 15 λ-windows × 250k steps forward and backward, both legs —
went from **~25–30 days to ~14–16 hours.** It now runs overnight.

---

## Is the speedup real, or is it cheating?

This is the question that matters for a scientific code. A faster integrator that
gives the wrong energy is worthless.

I did a single-point energy comparison on **identical coordinates** between the
old and new builds:

| term | old build | GPU-resident |
|------|-----------|--------------|
| ELECT | −226537.8616 | −226537.8756 |
| VDW | 19116.7494 | 19116.7546 |
| **POTENTIAL** | **−199582.4190** | **−199582.4277** |

Agreement to ~8 significant figures — a relative difference of 4×10⁻⁸, i.e.
pure float32 rounding. The speedup is **not** an accuracy tradeoff.

---

## Things that did *not* help

- **More GPU memory.** NAMD peaks at ~590 MiB on this system. The card is a 32 GB
  V100 SKU (16 GiB exposed to my container). A bigger-memory GPU would have changed
  **nothing** — the limit is compute/latency, not capacity.
- **Bigger `stepspercycle` / `pairlistsPerCycle` / `margin`.** No measurable effect.
- **More output throttling.** `restartfreq`/`dcdfreq`/`XSTFreq` at 500 measured
  identically to 5000.

---

## Practical gotchas worth knowing

1. **`GPUresident` in exactly one place.** Setting it both in the config *and* as
   `--GPUresident on` on the command line is fatal:
   `ERROR: Multiple definitions of 'GPUresident'`. Pick one mechanism per workflow.
2. **Don't misread `PERFORMANCE:` lines.** That number is a *cumulative* average from
   step 0. If your run starts with `minimize 5000`, it'll read ~3.9 ns/day during
   minimization and ~23 ns/day well into the MD phase, while the *instantaneous* rate
   (`TIMING:` sec/step) is 49.8. Judge throughput from `TIMING:`, not the average.
3. **`CUDASOAintegrate` is deprecated** — the modern spelling is `GPUresident on`.
4. **Lone pairs are disabled** under GPU-resident (`Disabling lonepair support due to
   incompatability with GPU-resident`). Harmless if your system has none — but
   re-check if you ever simulate a CGenFF halogen with lone pairs.

---

## The two-line takeaway

```bash
# 1. Build with single-node CUDA (not just --with-cuda):
./config Linux-x86_64-g++.gpuresident --charm-arch multicore-linux-x86_64 \
  --with-single-node-cuda --cuda-prefix /usr/local/cuda-11.8 \
  --cuda-gencode arch=compute_70,code=sm_70 && make -j32

# 2. Run with ONE PE, not +p8:
namd3 +p1 +devices 0 run.namd      # with `GPUresident on` in the config
```

*Hardware: NVIDIA Tesla V100-SXM2-32GB (sm_70), CUDA 11.8. System: 64,651-atom
protein–ligand complex. Software: NAMD 3.0.3, Charm++ multicore-linux-x86_64.*
