# 基于 NAMD 的相对结合自由能（RBFE）完整指南

> 内容包含三部分：① 为什么物理模拟在今天仍是药物发现的"金标准"；② Baidu AI Studio 上编译 NAMD 与每会话 `setup.sh` 的环境准备；③ 6I5I 双拓扑 FEP 的手动分步运行工作流。
> 约定：除代码/命令/路径外，正文为中文。代码块可直接复制执行。

---

## 1. 为什么今天还要跑基于物理的 RBFE？

一个绕不开的问题：AI 模型（深度学习打分函数、生成模型、AlphaFold 类方法）在药物发现里越来越强，为什么我们还要用 NAMD 花十几个小时算一个"甲基 → 氢"这样的小扰动？

**AI 快而广，物理模拟准而深 —— 两者定位完全不同。**

1. **AI 擅长"广撒网"，不擅长"精确裁决"。**
   AI 的本质是统计/关联模型，靠大量数据学规律，擅长在巨大的化学空间里快速生成候选、粗筛、给出"可能活性"的排序。但结合亲和力的实验数据稀少、噪声大、化学类型千差万别，ML 打分在"两个只差一个基团的类似物之间分辨 <1 kcal/mol 的差别"时往往不够可靠，对未见过的化学类型也容易失效（泛化差、甚至幻觉）。

2. **FEP 有严格的热力学根基，是"原理上精确"的裁决者。**
   自由能微扰（Zwanzig 公式 / BAR / MBAR）建立在统计力学之上：只要采样充分收敛、力场足够准确，结果在原理上就是精确的。它**不需要任何训练数据**，原则上可以迁移到任何能被参数化的分子。这正是药物化学里最难的决策场景——"去掉/加上一个基团，结合变好还是变坏、差多少"——AI 最不稳、而相对结合自由能（RBFE）最擅长的地方（良好系列上业界公认可达 ~1 kcal/mol 精度）。

3. **物理模拟给出"为什么"，AI 黑箱给不出。**
   FEP/MD 能把差异分解为范德华、静电、去溶剂化、构象熵等项，告诉你**是哪个相互作用在驱动差异**，从而指导下一步化学修饰。这是打分函数很难提供的信息。

4. **今天的范式是"AI 出主意，物理出结论"。**
   生成式 AI 负责扩大候选空间，FEP/MD 负责对接近的候选做**精确排序**，把需要合成与实验验证的分子压到最少。物理模拟依然是最后"拍板"所用的**金标准**；AI 的作用是让拍板的次数变少、变便宜。

回到本教程：要回答的不是"这个分子有没有活性"，而是"把 N–CH₃ 换成 N–H 之后结合自由能变化多少"——这是制药公司做**先导化合物优化（lead optimization）**时每天都在问的问题，也是 RBFE 真正的工业价值所在。

---

## 2. 使用环境：Baidu AI Studio

### 2.1 版本一览

| 组件 | 版本/位置 |
|------|-----------|
| NAMD | **3.0.3**（源码编译） |
| Charm++ | 8.0.0（打包在 NAMD 源码包里） |
| CUDA toolkit | 11.8（`/usr/local/cuda-11.8`，nvcc V11.8.89） |
| GPU / 架构 | Tesla V100-SXM2-32GB，sm_70 → `arch=compute_70,code=sm_70` |
| 编译目标 | `Linux-x86_64-g++.gpuresident` ← **生产用**（快 3.9 倍） |
| 生产二进制 | `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3` |

### 2.2 为什么需要一个 `setup.sh`（AI Studio 云环境的现实）

- Baidu AI Studio 是**临时性云端运行环境**：每次**重新启动项目都会得到一个全新的 shell/运行环境**，PATH、环境变量、系统软链接、会话内安装的软件包都会**丢失**。
- 只有**持久化目录 `/home/aistudio`** 里的文件能活下来（这也是为什么 NAMD 源码树和本仓库放在这里，但系统"记不住"它们）。
- NAMD 是**源码编译、非 `apt` 安装**，永远不会自动出现在 PATH 上；它的位置以及它依赖的 CUDA 根路径，每次都要**重新声明**。
- 云任务还会在启动约 **4.5 小时后被杀死**（见 `checkpoint_status.sh` 头部注释），于是你得反复重新启动计算会话——**每次重启都需要重新配好环境**。

所以用**一个 `setup.sh`**（每个会话一开始 `source` 一次）把 shell 重新指向 NAMD/CUDA，运行脚本就能原封不动地使用 `$NAMD`。

### 2.3 每个会话第一步：`source setup.sh`

```bash
# --- 任意目录
source /home/aistudio/work/setup.sh
```

它导出的变量（幂等，可反复执行）：

```bash
export NAMD_ROOT=/home/aistudio/NAMD_3.0.3_Source
export NAMD="$NAMD_ROOT/Linux-x86_64-g++.gpuresident/namd3"   # 生产二进制
export NAMD_CPU="$NAMD_ROOT/Linux-x86_64-g++/namd3"            # 旧版对照
export CUDA_HOME=/usr/local/cuda-11.8
export PATH="$NAMD_ROOT/Linux-x86_64-g++.gpuresident:$PATH"
```

source 后的自检：

```bash
$NAMD +p1 +devices 0   # 打印 NAMD 版本横幅即正常（无配置会随即退出）
echo "$NAMD"           # → /home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3
nvidia-smi             # 确认 GPU 空闲可用
```

若提示 `⚠ NAMD binary not found`，按 2.4 重新编译，或恢复源码树。

### 2.4 从源码编译 NAMD 3.0.3（一次性，约 15–25 分钟）

机器已具备：gcc/g++ 9.4、OpenMPI、CMake 3.24、conda Tcl、CUDA 11.8。**不需要 gfortran，不需要 sudo。**

#### ① 解压

```bash
cd /home/aistudio
tar xf NAMD_3.0.3_Source.tar.gz     # 内含 charm-8.0.0.tar
cd NAMD_3.0.3_Source
tar xf charm-8.0.0.tar
```

#### ② 编译 Charm++（multicore —— 单节点 GPU-resident 的正确架构）

```bash
cd charm-8.0.0
./build charm++ multicore-linux-x86_64 --with-production -j8
```

#### ③ 提供 Tcl 与 FFTW

```bash
# 从 http://www.ks.uiuc.edu/Research/namd/libraries/ 下载：
# tcl8.6.13-linux-x86_64.tar.gz、tcl8.6.13-linux-x86_64-threaded.tar.gz、fftw-linux-x86_64.tar.gz
```

⚠ **FFTW 必须用 `-fPIC` 重新编译。** 预编译的 `libsfftw.a`/`libsrfftw.a` 是 FFTW 2.1.x、**不带 `-fPIC`**，在 Ubuntu 20.04 的 PIE 默认下最终链接会失败：
`relocation R_X86_64_32 ... can not be used when making a PIE object`。

```bash
cd /home/aistudio
tar xf fftw-linux-x86_64.tar.gz
# 用 -fPIC 重建 FFTW 2.1.5，然后把产物拷进 NAMD_3.0.3_Source/fftw/：
#   libsfftw.a libsrfftw.a sfftw.h srfftw.h
# 并修改 srfftw.h：  #include "fftw.h"  →  #include <sfftw.h>
```

#### ④ 配置 + make —— GPU-resident（生产）

```bash
cd /home/aistudio/NAMD_3.0.3_Source
./config Linux-x86_64-g++.gpuresident \
  --charm-arch multicore-linux-x86_64 \
  --with-single-node-cuda \
  --cuda-prefix /usr/local/cuda-11.8 \
  --cuda-gencode arch=compute_70,code=sm_70
make -j32
```

`--with-single-node-cuda` 正是开启 `GPUresident` 的开关；只加 `--with-cuda` 无法跑 GPU-resident。
（如需旧版对照，把目标名换成 `Linux-x86_64-g++`、选项换成 `--with-cuda` 再编一份即可。）

#### ⑤ 验证

```bash
cd /home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident
./namd3 +p1 +devices 0 ../../src/alanin     # GPU 冒烟测试
```

### 2.5 运行规则（重要）

- **只用 GPU-resident 二进制，且 `+p1`（`+devices 0`）**。此 GPU 上 CPU 线程越多越慢（吞吐 ≈ 与 PE 数成反比）。
- **同时只跑一个 NAMD 任务**——只有一块 V100，绝不拆分。
- `GPUresident on` 写进配置文件，且**只能写在一处**（配置文件内，不要同时命令行传）。
- 每次 AI Studio 重启后 → `source ~/work/setup.sh` → 重跑相应生产腿（会自动从当前 λ-window 续跑）。

---

## 3. 6I5I 双拓扑 FEP —— 手动分步运行工作流

本节为**手动、逐条执行**的完整流程：每个代码块自上而下执行；两条相邻命令之间要**等上一条跑完**再执行。任何一步被中断都安全——重跑对应 `run` 命令会跳过已完成窗口、从进行中的窗口续跑。

### 3.1 运行状态

本示例**已跑完**：complex / solvent 双腿的 forward + backward 全部 15 个窗口
（每窗口 250 000 步 / 500 ps）都已完成并组装。最终结果（见 §3.8）：

| 量 | 值 |
|----|----|
| dG_complex | +3.882 kcal/mol |
| dG_solvent | +3.988 kcal/mol |
| **ΔΔG** | **−0.106 kcal/mol** |
| 95% 置信区间（moving-block bootstrap） | [−0.355, +0.118] |

即**与 0 无法区分**——15 × 500 ps 的采样量不足以分辨这对配体，
详见根目录 [README §5](README.md#5-worked-example--6i5i-clk1) 的说明。

> 实时状态用下面命令查询，别依赖上面的静态表：
> ```bash
> cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP
> python3 fep_run.py windows complex md_backward   # 返回 JSON 状态
> ```

### 3.2 常量

```bash
source /home/aistudio/work/setup.sh   # 每次 AI Studio 会话先执行：导出 $NAMD、$CUDA_HOME、PATH
ROOT=/home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP
NAMD=${NAMD:-/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3}
```

- 二进制是 **GPU-resident** 的 namd3，永远用 `+p1 +devices 0`。
- 环境与编译背景见第 2 节。

### 3.3 目录结构

```
6I5I_DUAL_FEP/
├── prepare_hybrid.py        # 配体映射 → hybrid.rtf/prm/pdb        （§3.4 可选重建）
├── build_system.py          # psfgen → solvate → ionize            （§3.4，需 VMD）
├── write_fep_inputs.py      # ionized.fep B-factor + *.namd + fep.tcl（§3.4）
├── fep_run.py               # 逐 λ-window NAMD 驱动（核心工具）
│   #   fep_run.py run      <leg> <forward|backward>   [--start N] [--from STEM] [--steps N]
│   #   fep_run.py assemble <leg> <md_forward|md_backward>
│   #   fep_run.py windows  <leg> <md_forward|md_backward>      # 状态探测 → JSON
├── complex/                 # 64,651 原子：蛋白 + 配体 + 水
│   ├── nvt_equil.namd  npt_equil.namd
│   ├── md_forward.namd  md_backward.namd      # 单进程模板（配置来源）
│   ├── md_{forward,backward}_wNN.namd         # 逐窗配置（由 fep_run 生成）
│   ├── md_{forward,backward}_wNN.{coor,vel,xsc,fepout,dcd,log}
│   ├── .md_*_wNN.done                          # 完成标记（跳过逻辑依据）
│   ├── md_forward_combined.fepout              # assemble 输出
│   └── window_snapshots/
├── solvent/                 # 5,013 原子：仅配体 + 水（文件布局同 complex）
├── run_all.sh               # ⚠ 已废弃的单进程流水线 —— 不要再运行
└── fep.tcl                  # 16 个 λ 的调度表（与 fep_run.py 内 LAM[] 一致）
```

### 3.4 Phase 1 —（可选）重新生成输入，当前已完成

```bash
# --- from: $ROOT
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP
python3 prepare_hybrid.py      # 配体映射 -> hybrid/
python3 build_system.py        # complex/ + solvent/（需 VMD：psfgen、solvate、ionize）
python3 write_fep_inputs.py    # ionized.fep B-factor、fep.tcl、*.namd
```

⚠ 重跑会**覆盖**当前已平衡好的体系，只在想"推倒重来"时执行。

### 3.5 Phase 2 — 平衡（每个体系：先 nvt 再 npt）

nvt 负责赋初速度，npt 从 nvt 续跑。complex、solvent 各做一遍。

```bash
# --- from: $ROOT/complex
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/complex
$NAMD +p1 +devices 0 nvt_equil.namd > nvt_equil.log 2>&1 && touch .nvt_equil.done
$NAMD +p1 +devices 0 npt_equil.namd > npt_equil.log 2>&1 && touch .npt_equil.done
```

```bash
# --- from: $ROOT/solvent
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/solvent
$NAMD +p1 +devices 0 nvt_equil.namd > nvt_equil.log 2>&1 && touch .nvt_equil.done
$NAMD +p1 +devices 0 npt_equil.namd > npt_equil.log 2>&1 && touch .npt_equil.done
```

提醒：这类从重启文件续跑的配置**不能写 `temperature`**（npt/md 靠 `binvelocities` 继承初速度）。
每个方向第 0 个窗口以 `npt_equil` 为起点。

### 3.6 Phase 3 — 生产 FEP（核心）

每个 λ-window 单独启动一次 NAMD、窗口间首尾相接（后窗从前窗的 coor/vel/xsc 续跑）。
顺序固定为：**complex forward → complex backward → solvent forward → solvent backward**。
complex 每窗约 17 分钟，solvent 每窗约 12 分钟。

```bash
# --- from: $ROOT（fep_run.py 需在此目录运行）
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP

python3 fep_run.py run complex forward     # ← 已完成，会全部 skip 后返回
python3 fep_run.py run complex backward    # ← 当前续跑点：重跑进行中的 w10，再跑 w11..14
python3 fep_run.py run solvent forward     # ← 15 个窗口
python3 fep_run.py run solvent backward    # ← 15 个窗口
```

说明：
- 已完成的窗口靠 `.md_*_wNN.done` 标记自动跳过——可随时重跑。
- 被打断的窗口会从**上一窗口边界**重做（最多损失 ≤ ~18 分钟）。
- 可选自定义起点：
  `python3 fep_run.py run complex backward --start 10`
  `python3 fep_run.py run complex forward --start 6 --from window_snapshots/md_forward_w06`

随时探测状态（JSON）：

```bash
python3 fep_run.py windows complex md_forward
python3 fep_run.py windows complex md_backward
python3 fep_run.py windows solvent md_forward
python3 fep_run.py windows solvent md_backward
```

### 3.7 Phase 4 — 单窗口裸跑 NAMD（可选 / 诊断用）

`fep_run.py run` 内部就是下面这条命令、一个窗口一次。只有当你手上有现成的逐窗 `.namd` 时才直接跑（`.namd` 是 `run` 执行时现场生成的）：

```bash
# --- from: $ROOT/complex（或 solvent/）
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/complex
$NAMD +p1 +devices 0 md_backward_w10.namd > md_backward_w10.log 2>&1 && touch .md_backward_w10.done
```

第 0 个窗口以 `npt_equil` 为续跑起点，第 N 个窗口从 `md_*_w{N-1}` 续跑。

### 3.8 Phase 5 — 组装 + BAR 分析

#### 组装每个腿的 combined fepout（在其 `run` 结束后）

```bash
# --- from: $ROOT
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP
python3 fep_run.py assemble complex md_forward     # → complex/md_forward_combined.fepout
python3 fep_run.py assemble complex md_backward
python3 fep_run.py assemble solvent md_forward
python3 fep_run.py assemble solvent md_backward
```

#### BAR 求 ΔΔG

```bash
# --- from: /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP

# 推荐：完整审计（解析器校验 + BAR + 迟滞 + 平稳性 + block bootstrap 误差棒）
python3 audit_fep.py

# 只要一个数：
python3 analyze_fep.py bar \
   complex/md_forward_combined.fepout complex/md_backward_combined.fepout \
   solvent/md_forward_combined.fepout solvent/md_backward_combined.fepout
```

> ⚠️ 分析**只用本仓库的 `analyze_fep.py` / `audit_fep.py`（推荐后者）**，
> 它们直接解析 NAMD 的 `FepEnergy` 原始行。任何靠正则去抓
> "Free energy change" 文本的做法，抓到的都是累计值 `net change until now`
> 而不是每窗口的值：实测给出 −16.821（真值 −5.455），而且**看起来像个正常结果**。
>
> ⚠️ `analyze_fep.py` / `audit_fep.py` 默认会**丢弃每个窗口的前 99 个样本**
> （`alchEquilSteps 50000` ÷ `alchOutFreq 500` 的平衡段）。NAMD 自身的
> `dE_avg`/`dG` 列也只统计平衡之后的 401 个样本。想复现旧行为用 `--no-trim`。

### 3.9 监控与恢复规则

```bash
# 实时步数输出（指向正在跑的窗口）：
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP && tail -f complex/md_backward_w10.log
# 每 2 分钟一行状态快照（只读监控，可选）：
cd /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP && setsid nohup bash checkpoint_status.sh 120 >/dev/null 2>&1 < /dev/null &
```

- **速率看 `TIMING:` 行**：complex 约 41 ns/day，solvent 约 59 ns/day（忽略 `PERFORMANCE:` 的累计均值）。
- **被杀/续跑**：本机任务约在启动 4.5 小时后被杀。任何中断后，直接重跑第 3.6 节对应的 `run` 行即可续跑。
- **健康检查**：`nvidia-smi` 应显示 ~700 MiB 且利用率高；日志尾部持续出现 `TIMING`/`ENERGY`/`FEP:` 与周期性的 `WRITING ... RESTART`。

---

## 4. 结束语

把"AI 生成候选 + 物理模拟裁决"两条腿都跑通，你就拥有了和大型制药公司 FEP 平台一致的决策流程：**用统计力学回答"改一下结合差多少"，用最少的合成实验拿到最可靠的答案。** 祝跑批顺利！
