"""Tests for the upload contract.

Run:  python3 -m unittest discover -s tests -v
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from fep_web import contract          # noqa: E402


def _sdf(nheavy: int, nhydro: int, title: str = "lig") -> str:
    n = nheavy + nhydro
    out = [title, "  test", "",
           f"{n:3d}{0:3d}  0  0  0  0  0  0  0  0999 V2000"]
    for _ in range(nheavy):
        out.append(f"{0.0:10.4f}{0.0:10.4f}{0.0:10.4f} C   0  0  0  0  0  0  0"
                   f"  0  0  0  0  0  0")
    for _ in range(nhydro):
        out.append(f"{0.0:10.4f}{0.0:10.4f}{0.0:10.4f} H   0  0  0  0  0  0  0"
                   f"  0  0  0  0  0  0")
    out += ["M  END", "$$$$"]
    return "\n".join(out) + "\n"


def _pdb(nheavy: int, nhydro: int, chains=("A",)) -> str:
    out = []
    serial = 0
    for ci, ch in enumerate(chains):
        for i in range(nheavy // len(chains)):
            serial += 1
            out.append(f"ATOM  {serial:5d}  CA  ALA {ch}{i + 1:4d}    "
                       f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
        for i in range(nhydro // len(chains)):
            serial += 1
            out.append(f"ATOM  {serial:5d}  HA  ALA {ch}{i + 1:4d}    "
                       f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           H")
    out.append("END")
    return "\n".join(out) + "\n"


def _system(tmp: Path, *, manifest=None, protein=None, ref=None, mut=None,
            name="SYS") -> Path:
    root = tmp / name
    root.mkdir(parents=True, exist_ok=True)
    (root / "system.json").write_text(
        json.dumps({} if manifest is None else manifest))
    (root / "protein.pdb").write_text(
        _pdb(200, 150) if protein is None else protein)
    (root / "ref.sdf").write_text(_sdf(20, 12) if ref is None else ref)
    (root / "mut.sdf").write_text(_sdf(19, 11) if mut is None else mut)
    return root


class TestValidInput(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_minimal_valid_system_passes(self):
        v = contract.validate(_system(self.tmp))
        self.assertTrue(v.ok, v.summary())
        self.assertEqual(v.errors, [])
        self.assertEqual(set(v.ligands), {"ref", "mut"})

    def test_defaults_are_filled_in(self):
        v = contract.validate(_system(self.tmp))
        self.assertEqual(v.manifest["segments"], "auto")
        self.assertEqual(v.manifest["padding"], 15.0)
        self.assertEqual(v.manifest["temperature"], 300.0)

    def test_explicit_manifest_keys_win(self):
        root = _system(self.tmp, manifest={"segments": [["P1", 1, 400]],
                                           "temperature": 310.0})
        v = contract.validate(root)
        self.assertTrue(v.ok, v.summary())
        self.assertEqual(v.manifest["segments"], [["P1", 1, 400]])
        self.assertEqual(v.manifest["temperature"], 310.0)

    def test_unknown_manifest_key_is_a_warning_not_an_error(self):
        v = contract.validate(_system(self.tmp, manifest={"wibble": 1}))
        self.assertTrue(v.ok)
        self.assertTrue(any("wibble" in str(p) for p in v.warnings))


class TestRejections(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_missing_manifest(self):
        root = _system(self.tmp)
        (root / "system.json").unlink()
        v = contract.validate(root)
        self.assertFalse(v.ok)
        self.assertTrue(any(contract.MANIFEST in str(p) for p in v.errors))

    def test_invalid_json_manifest(self):
        root = _system(self.tmp)
        (root / "system.json").write_text("{not json")
        v = contract.validate(root)
        self.assertFalse(v.ok)
        self.assertIn("invalid JSON", v.summary())

    def test_missing_protein(self):
        root = _system(self.tmp)
        (root / "protein.pdb").unlink()
        self.assertFalse(contract.validate(root).ok)

    def test_protein_without_hydrogens_is_allowed_with_a_warning(self):
        """psfgen adds them from the topology.

        The 6I5I protein.pdb that produced the published result has zero
        hydrogens, so rejecting it would reject the pipeline's own working
        input. Protonation state is decided by build_system.py's pdbalias.
        """
        v = contract.validate(_system(self.tmp, protein=_pdb(200, 0)))
        self.assertTrue(v.ok, v.summary())
        self.assertTrue(any("psfgen will add them" in str(p) for p in v.warnings))
        self.assertEqual(v.errors, [])

    def test_tiny_protein_is_rejected(self):
        v = contract.validate(_system(self.tmp, protein=_pdb(5, 3)))
        self.assertFalse(v.ok)

    def test_missing_ligand(self):
        root = _system(self.tmp)
        (root / "mut.sdf").unlink()
        v = contract.validate(root)
        self.assertFalse(v.ok)
        self.assertIn("mut", v.summary())

    def test_ligand_without_hydrogens_is_rejected(self):
        """The silent killer: heavy-atom-only input gives wrong parameters."""
        v = contract.validate(_system(self.tmp, ref=_sdf(20, 0)))
        self.assertFalse(v.ok)
        self.assertIn("no hydrogens", v.summary())

    def test_truncated_ligand_is_rejected(self):
        v = contract.validate(_system(self.tmp, ref=_sdf(3, 2)))
        self.assertFalse(v.ok)
        self.assertIn("truncated", v.summary())

    def test_unparseable_ligand_is_rejected(self):
        v = contract.validate(_system(self.tmp, ref="this is not a molecule\n"))
        self.assertFalse(v.ok)
        self.assertIn("cannot parse", v.summary())

    def test_missing_directory(self):
        self.assertFalse(contract.validate(self.tmp / "nope").ok)


class TestWarnings(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_identical_ligands_warn(self):
        v = contract.validate(_system(self.tmp, ref=_sdf(20, 12),
                                      mut=_sdf(20, 12)))
        self.assertTrue(v.ok)
        self.assertTrue(any("identical atom counts" in str(p) for p in v.warnings))

    def test_multichain_protein_warns_about_segments(self):
        v = contract.validate(_system(self.tmp, protein=_pdb(200, 150,
                                                             chains=("A", "B"))))
        self.assertTrue(v.ok)
        self.assertTrue(any("chains detected" in str(p) for p in v.warnings))


class TestFormats(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_mol2_ligand_accepted(self):
        root = _system(self.tmp)
        (root / "ref.sdf").unlink()
        atoms = ["@<TRIPOS>MOLECULE", "lig", " 10 0 0 0 0", "",
                 "@<TRIPOS>ATOM"]
        atoms += [f"{i:7d} C{i:<4d} 0.0 0.0 0.0 C.3  1  LIG  0.0"
                  for i in range(1, 8)]
        atoms += [f"{i:7d} H{i:<4d} 0.0 0.0 0.0 H    1  LIG  0.0"
                  for i in range(8, 11)]
        (root / "ref.mol2").write_text("\n".join(atoms) + "\n")
        v = contract.validate(root)
        self.assertTrue(v.ok, v.summary())
        self.assertEqual(v.ligands["ref"].suffix, ".mol2")

    def test_sdf_counts_line_is_parsed(self):
        p = self.tmp / "x.sdf"
        p.write_text(_sdf(20, 12))
        self.assertEqual(contract.read_ligand(p), (20, 12))

    def test_mol2_counts(self):
        p = self.tmp / "x.mol2"
        p.write_text("@<TRIPOS>MOLECULE\nm\n\n@<TRIPOS>ATOM\n"
                     "1 C1 0 0 0 C.3 1 L 0\n2 H1 0 0 0 H 1 L 0\n")
        self.assertEqual(contract.read_ligand(p), (1, 1))

    def test_bad_sdf_is_rejected_by_reader(self):
        p = self.tmp / "x.sdf"
        p.write_text("nope\n")
        with self.assertRaises(ValueError):
            contract.read_ligand(p)

    def test_unsupported_extension(self):
        with self.assertRaises(ValueError):
            contract.read_ligand(self.tmp / "x.xyz")


class TestAgainstRealData(unittest.TestCase):
    """Validate the contract against the repo's own real ligand files."""

    LIG = REPO / "6I5I_FEP" / "1_ligand"

    def test_real_6i5i_ligands_have_hydrogens(self):
        ref = self.LIG / "ref_crystal.sdf"
        mut = self.LIG / "mut_aligned.sdf"
        if not ref.exists() or not mut.exists():
            self.skipTest("6I5I_FEP ligand files not present")
        ref_h, ref_hyd = contract.read_ligand(ref)
        mut_h, mut_hyd = contract.read_ligand(mut)
        self.assertGreater(ref_hyd, 0, "ref_crystal.sdf has no explicit H")
        self.assertGreater(mut_hyd, 0, "mut_aligned.sdf has no explicit H")
        self.assertGreater(ref_h, mut_h, "desmethyl should have fewer heavy atoms")


if __name__ == "__main__":
    unittest.main(verbosity=2)
