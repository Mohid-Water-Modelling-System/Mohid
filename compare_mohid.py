#!/usr/bin/env python3
"""
compare_mohid.py — MOHID output validation tool
================================================
Compares HDF5 results and text output files (*.dat, *.srr) between a baseline
and a new run, reporting maximum differences and NaN/Inf anomalies.

USAGE — auto mode (detects *_original / *_Original pairs in the res folder):
    python compare_mohid.py D:\\...\\res

USAGE — explicit mode (compare two files or two directories by matching name):
    python compare_mohid.py D:\\...\\res\\RunOff_55_original.hdf5  D:\\...\\res\\RunOff_55.hdf5
    python compare_mohid.py D:\\...\\res\\Run55_original           D:\\...\\res\\Run55

OPTIONS:
    --tol-abs FLOAT   Absolute tolerance for FAIL  (default: 1e-5)
    --tol-rel FLOAT   Relative tolerance for FAIL  (default: 1e-4)
                      A dataset fails if BOTH abs AND rel are exceeded.
    --verbose, -v     Show the top differing datasets / file sections
    --no-color        Disable ANSI colours

REQUIREMENTS:
    pip install h5py numpy
"""

import argparse
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(f"Missing dependency: {exc}\nInstall with:  pip install h5py numpy")
    sys.exit(2)

# Force UTF-8 output so box-drawing chars and ✓/✗ work in terminals and
# file redirection on Windows (Python 3.7+).
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ── colour helpers ────────────────────────────────────────────────────────────

def _color_ok() -> bool:
    if "--no-color" in sys.argv:
        return False
    if sys.stdout.isatty():
        # Enable ANSI on Windows consoles that support it
        if os.name == "nt":
            os.system("")          # triggers VT100 mode
        return True
    return False

_USE_COLOR = _color_ok()

def _c(text: str, code: str) -> str:
    return f"\033[{code}m{text}\033[0m" if _USE_COLOR else text

def green(t):  return _c(t, "32")
def red(t):    return _c(t, "31")
def yellow(t): return _c(t, "33")
def bold(t):   return _c(t, "1")
def cyan(t):   return _c(t, "36")

SEP  = "-" * 80
SEP2 = "=" * 80

def fmt(v: float) -> str:
    """Format a float for the report (zero shows as 0, others as sci notation)."""
    if v == 0.0:
        return "0"
    return f"{v:.3e}"

# ── data classes ──────────────────────────────────────────────────────────────

@dataclass
class DatasetDiff:
    """Result for one HDF5 dataset or one named section of a text file."""
    name: str
    max_abs: float = 0.0
    max_rel: float = 0.0
    nan_new: int   = 0          # NaN/Inf counts in new file
    nan_base: int  = 0          # NaN/Inf counts in baseline
    status: str    = "pass"     # pass | fail | skip | shape_mismatch | missing
    note: str      = ""


@dataclass
class FileResult:
    label: str          # display name, e.g. "RunOff_55.hdf5"
    baseline: Path
    new: Path
    file_type: str      # hdf5 | text | unknown
    status: str = "pass"   # pass | fail | error | missing_new | missing_base
    diffs: List[DatasetDiff] = field(default_factory=list)
    error: str = ""

    @property
    def max_abs(self) -> float:
        return max((d.max_abs for d in self.diffs
                     if d.status not in ("skip", "missing", "shape_mismatch")), default=0.0)

    @property
    def max_rel(self) -> float:
        return max((d.max_rel for d in self.diffs
                     if d.status not in ("skip", "missing", "shape_mismatch")), default=0.0)

    @property
    def nan_new_total(self) -> int:
        return sum(d.nan_new for d in self.diffs)

    @property
    def nan_base_total(self) -> int:
        return sum(d.nan_base for d in self.diffs)

    @property
    def n_compared(self) -> int:
        return sum(1 for d in self.diffs if d.status not in ("skip",))


# ── HDF5 comparison ───────────────────────────────────────────────────────────

def _compare_hdf5(baseline: Path, new: Path, tol_abs: float, tol_rel: float) -> FileResult:
    result = FileResult(
        label=new.name, baseline=baseline, new=new, file_type="hdf5"
    )
    try:
        with h5py.File(baseline, "r") as fb, h5py.File(new, "r") as fn:
            # Collect all dataset paths from both files
            base_ds: dict[str, tuple] = {}
            new_ds:  set[str]         = set()

            def _collect_base(name, obj):
                if isinstance(obj, h5py.Dataset):
                    base_ds[name] = obj.shape

            def _collect_new(name, obj):
                if isinstance(obj, h5py.Dataset):
                    new_ds.add(name)

            fb.visititems(_collect_base)
            fn.visititems(_collect_new)

            for ds_path in sorted(base_ds):
                dd = _compare_hdf5_dataset(fb, fn, ds_path, tol_abs, tol_rel, new_ds)
                result.diffs.append(dd)
                if dd.status in ("fail", "missing", "shape_mismatch"):
                    result.status = "fail"

            # datasets present in new but not in baseline
            for ds_path in sorted(new_ds - set(base_ds)):
                result.diffs.append(
                    DatasetDiff(name=ds_path, status="missing",
                                note="present in new, absent in baseline")
                )

    except Exception as exc:
        result.status = "error"
        result.error  = str(exc)

    return result


def _compare_hdf5_dataset(
    fb: h5py.File, fn: h5py.File,
    path: str,
    tol_abs: float, tol_rel: float,
    new_ds: set,
) -> DatasetDiff:
    dd = DatasetDiff(name=path)

    if path not in new_ds:
        dd.status = "missing"
        dd.note   = "present in baseline, absent in new"
        return dd

    ds_b = fb[path]
    ds_n = fn[path]

    # Skip non-numeric datasets
    if ds_b.dtype.kind in ("S", "U", "O"):
        dd.status = "skip"
        dd.note   = "string dataset"
        return dd

    if ds_b.shape != ds_n.shape:
        dd.status = "shape_mismatch"
        dd.note   = f"baseline{ds_b.shape} vs new{ds_n.shape}"
        return dd

    if ds_b.size == 0:
        dd.status = "skip"
        dd.note   = "empty"
        return dd

    a = ds_b[()].astype(np.float64)
    b = ds_n[()].astype(np.float64)

    dd.nan_base = int(np.sum(~np.isfinite(a)))
    dd.nan_new  = int(np.sum(~np.isfinite(b)))

    mask = np.isfinite(a) & np.isfinite(b)
    if np.any(mask):
        diff  = np.abs(a[mask] - b[mask])
        dd.max_abs = float(np.max(diff))
        dd.max_rel = float(np.max(diff / np.maximum(np.abs(a[mask]), 1e-300)))

    new_nan_appeared = bool(np.any(np.isfinite(a) & ~np.isfinite(b)))
    over_tol         = dd.max_abs > tol_abs and dd.max_rel > tol_rel
    if new_nan_appeared or over_tol:
        dd.status = "fail"

    return dd


# ── Text file comparison (*.dat, *.srr) ───────────────────────────────────────

_FLOAT_RE = re.compile(
    r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?"
)

def _extract_floats(path: Path) -> Tuple[np.ndarray, int]:
    """Return all numeric values in a text file as a float64 array + line count."""
    values  = []
    n_lines = 0
    with open(path, "r", encoding="latin-1", errors="replace") as fh:
        for line in fh:
            n_lines += 1
            for tok in _FLOAT_RE.findall(line):
                try:
                    values.append(float(tok.replace("D", "E").replace("d", "e")))
                except ValueError:
                    pass
    return np.asarray(values, dtype=np.float64), n_lines


def _compare_text(baseline: Path, new: Path, tol_abs: float, tol_rel: float) -> FileResult:
    result = FileResult(
        label=new.name, baseline=baseline, new=new, file_type="text"
    )
    try:
        arr_b, lines_b = _extract_floats(baseline)
        arr_n, lines_n = _extract_floats(new)

        if arr_b.size != arr_n.size:
            dd = DatasetDiff(
                name="(all values)",
                status="fail",
                note=(f"different value count — baseline: {arr_b.size:,}  "
                      f"new: {arr_n.size:,}  "
                      f"(baseline lines: {lines_b}, new lines: {lines_n})"),
            )
            result.diffs.append(dd)
            result.status = "fail"
            return result

        dd = DatasetDiff(name="(all values)")
        dd.nan_base = int(np.sum(~np.isfinite(arr_b)))
        dd.nan_new  = int(np.sum(~np.isfinite(arr_n)))

        mask = np.isfinite(arr_b) & np.isfinite(arr_n)
        if np.any(mask):
            diff = np.abs(arr_b[mask] - arr_n[mask])
            dd.max_abs = float(np.max(diff))
            dd.max_rel = float(np.max(diff / np.maximum(np.abs(arr_b[mask]), 1e-300)))

        new_nan_appeared = bool(np.any(np.isfinite(arr_b) & ~np.isfinite(arr_n)))
        over_tol         = dd.max_abs > tol_abs and dd.max_rel > tol_rel
        if new_nan_appeared or over_tol:
            dd.status = "fail"
            result.status = "fail"

        result.diffs.append(dd)

    except Exception as exc:
        result.status = "error"
        result.error  = str(exc)

    return result


# ── Auto-discovery within a single res directory ──────────────────────────────

# Suffixes that mark a saved baseline (case-insensitive)
_BASELINE_SUFFIXES = ("_original", "_phase1", "_baseline")

def _strip_baseline_suffix(name: str) -> Optional[str]:
    """If name ends with a known baseline suffix, return the stem without it."""
    lower = name.lower()
    for suf in _BASELINE_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return None


def _auto_discover(res_dir: Path) -> List[Tuple[Path, Path, str]]:
    """
    Return list of (baseline, new, label) pairs found in res_dir.
    Looks for:
      - Files: <stem>_original.ext   →   <stem>.ext
      - Dirs:  <name>_original       →   <name>
    """
    pairs: List[Tuple[Path, Path, str]] = []

    for item in sorted(res_dir.iterdir()):
        if item.is_dir():
            stem = _strip_baseline_suffix(item.name)
            if stem is None:
                continue
            new_path = res_dir / stem
            if new_path.is_dir():
                pairs.append((item, new_path, f"{stem}/ (directory)"))
        elif item.is_file():
            # e.g.  RunOff_55_original.hdf5  →  RunOff_55.hdf5
            suffix    = item.suffix                          # .hdf5
            bare_stem = item.stem                            # RunOff_55_original
            clean_stem = _strip_baseline_suffix(bare_stem)
            if clean_stem is None:
                continue
            new_path = res_dir / (clean_stem + suffix)
            if new_path.is_file():
                pairs.append((item, new_path, clean_stem + suffix))

    return pairs


# ── Compare a matched pair (file or directory) ────────────────────────────────

_HDF_EXT  = {".hdf5", ".h5", ".hdf"}
_TEXT_EXT = {".dat", ".srr"}


def _compare_pair(
    baseline: Path, new: Path,
    tol_abs: float, tol_rel: float,
    label: str,
) -> List[FileResult]:
    """Compare one matched pair; returns one or more FileResult objects."""

    if baseline.is_file() and new.is_file():
        return [_compare_single_file(baseline, new, tol_abs, tol_rel)]

    if baseline.is_dir() and new.is_dir():
        return _compare_directories(baseline, new, tol_abs, tol_rel)

    # Mismatched types
    fr = FileResult(
        label=label, baseline=baseline, new=new, file_type="unknown",
        status="error", error="one path is a file, the other a directory",
    )
    return [fr]


def _compare_single_file(
    baseline: Path, new: Path, tol_abs: float, tol_rel: float
) -> FileResult:
    ext = baseline.suffix.lower()
    if ext in _HDF_EXT:
        return _compare_hdf5(baseline, new, tol_abs, tol_rel)
    if ext in _TEXT_EXT:
        return _compare_text(baseline, new, tol_abs, tol_rel)
    fr = FileResult(
        label=new.name, baseline=baseline, new=new,
        file_type="unknown", status="pass",
    )
    fr.diffs.append(DatasetDiff(name="(file)", status="skip",
                                note=f"unsupported extension '{ext}'"))
    return fr


def _compare_directories(
    baseline_dir: Path, new_dir: Path, tol_abs: float, tol_rel: float
) -> List[FileResult]:
    results: List[FileResult] = []

    base_files = {
        f.name: f
        for f in baseline_dir.iterdir()
        if f.is_file() and not f.name.lower().startswith("stormwater")
    }
    new_files  = {
        f.name: f
        for f in new_dir.iterdir()
        if f.is_file() and not f.name.lower().startswith("stormwater")
    }
    all_names  = sorted(base_files.keys() | new_files.keys())

    for name in all_names:
        ext = Path(name).suffix.lower()
        if ext not in (_HDF_EXT | _TEXT_EXT):
            continue   # skip .fin, .map, .ver etc.

        label = f"{new_dir.name}/{name}"
        if name not in base_files:
            results.append(FileResult(
                label=label, baseline=baseline_dir / name,
                new=new_files[name], file_type="unknown",
                status="missing_base",
            ))
        elif name not in new_files:
            results.append(FileResult(
                label=label, baseline=base_files[name],
                new=new_dir / name, file_type="unknown",
                status="missing_new",
            ))
        else:
            fr = _compare_single_file(base_files[name], new_files[name], tol_abs, tol_rel)
            fr.label = label
            results.append(fr)

    return results


# ── Reporting ─────────────────────────────────────────────────────────────────

def _status_badge(status: str) -> str:
    return {
        "pass":         green("PASS [OK]"),
        "fail":         red("FAIL [!!]"),
        "error":        red("ERROR"),
        "missing_new":  yellow("WARN - missing in new"),
        "missing_base": yellow("WARN - missing in baseline"),
    }.get(status, status)


def _print_file_result(fr: FileResult, verbose: bool, tol_abs: float, tol_rel: float) -> None:
    print(f"\n{SEP}")
    print(f"{bold(fr.label)}   [{_status_badge(fr.status)}]")
    print(f"  baseline : {fr.baseline}")
    print(f"  new      : {fr.new}")

    if fr.status in ("missing_new", "missing_base"):
        return
    if fr.status == "error":
        print(f"  {red('error')} : {fr.error}")
        return

    print(f"  type     : {fr.file_type}")
    print(f"  compared : {fr.n_compared} dataset(s)/section(s)")
    print(f"  max abs  : {fmt(fr.max_abs)}")
    print(f"  max rel  : {fmt(fr.max_rel)}")

    if fr.nan_base_total:
        print(f"  NaN/Inf baseline : {fr.nan_base_total}")
    if fr.nan_new_total:
        nan_str = red(str(fr.nan_new_total)) if fr.nan_new_total > fr.nan_base_total else str(fr.nan_new_total)
        print(f"  NaN/Inf new      : {nan_str}")

    if verbose or fr.status == "fail":
        failing = [d for d in fr.diffs if d.status == "fail"]
        shape_issues = [d for d in fr.diffs if d.status == "shape_mismatch"]
        missing = [d for d in fr.diffs if d.status == "missing"]
        notable = failing + shape_issues + missing

        if notable:
            print(f"\n  {bold('Differences')} (abs > {tol_abs:.0e} AND rel > {tol_rel:.0e}):")
            for d in notable[:30]:
                line = f"    {d.name}"
                if d.status == "fail":
                    line += f"  abs={fmt(d.max_abs)}  rel={fmt(d.max_rel)}"
                    if d.nan_new > d.nan_base:
                        line += f"  {red(f'NaN/Inf+{d.nan_new - d.nan_base}')}"
                else:
                    line += f"  [{d.status}] {d.note}"
                print(line)
            if len(notable) > 30:
                print(f"    … and {len(notable) - 30} more")

        # Show top-5 worst passing datasets if verbose
        if verbose and fr.file_type == "hdf5":
            passing = sorted(
                [d for d in fr.diffs if d.status == "pass" and d.max_abs > 0],
                key=lambda d: d.max_abs, reverse=True,
            )[:5]
            if passing:
                print(f"\n  {bold('Top diffs (passing datasets):')} ")
                for d in passing:
                    print(f"    {d.name}  abs={fmt(d.max_abs)}  rel={fmt(d.max_rel)}")


def _print_summary(all_results: List[FileResult]) -> None:
    print(f"\n{SEP2}")
    n_total = len(all_results)
    n_pass  = sum(1 for r in all_results if r.status == "pass")
    n_fail  = sum(1 for r in all_results if r.status == "fail")
    n_warn  = sum(1 for r in all_results if r.status in ("missing_new", "missing_base"))
    n_err   = sum(1 for r in all_results if r.status == "error")

    if n_fail == 0 and n_err == 0 and n_warn == 0:
        verdict = green(f"ALL {n_pass}/{n_total} FILES PASS — outputs are equivalent")
    else:
        parts = [f"{n_pass}/{n_total} pass"]
        if n_fail: parts.append(red(f"{n_fail} FAIL"))
        if n_warn: parts.append(yellow(f"{n_warn} warning"))
        if n_err:  parts.append(red(f"{n_err} error"))
        verdict = "   ".join(parts)

    print(f"SUMMARY  {verdict}")
    print(SEP2)


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare MOHID output files (HDF5, .dat, .srr) between a baseline and a new run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (auto mode — detects *_original / *_phase1 pairs):
  python compare_mohid.py D:\\...\\Performance_LargeModel_quickversion\\res
  python compare_mohid.py D:\\...\\Performance_LargeModel_quickversion_noRain\\res

Examples (explicit mode):
  python compare_mohid.py D:\\...\\res\\RunOff_55_original.hdf5  D:\\...\\res\\RunOff_55.hdf5
  python compare_mohid.py D:\\...\\res\\Run55_original           D:\\...\\res\\Run55
  python compare_mohid.py D:\\...\\res\\Run55_phase1             D:\\...\\res\\Run55
""",
    )
    parser.add_argument("path_a", help="Baseline file/directory  OR  res directory (auto mode)")
    parser.add_argument("path_b", nargs="?", default=None,
                        help="New file/directory (omit for auto mode)")
    parser.add_argument("--tol-abs", type=float, default=1e-5,
                        help="Absolute tolerance (default: 1e-5). Both abs AND rel must exceed to FAIL.")
    parser.add_argument("--tol-rel", type=float, default=1e-4,
                        help="Relative tolerance (default: 1e-4).")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show per-dataset detail and top-5 worst diffs")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable ANSI colour output")
    args = parser.parse_args()

    if args.no_color:
        global _USE_COLOR
        _USE_COLOR = False

    path_a = Path(args.path_a)
    if not path_a.exists():
        print(f"Error: path not found: {path_a}", file=sys.stderr)
        sys.exit(2)

    # ── Determine mode ────────────────────────────────────────────────────────

    auto_mode = (args.path_b is None)

    if auto_mode:
        # path_a must be a directory (the res folder)
        if not path_a.is_dir():
            print("Error: in auto mode the argument must be a directory.", file=sys.stderr)
            sys.exit(2)
        pairs = _auto_discover(path_a)
        if not pairs:
            print(
                f"No baseline pairs found in: {path_a}\n"
                f"Expected files/dirs with suffix: {_BASELINE_SUFFIXES}",
                file=sys.stderr,
            )
            sys.exit(2)
        print(SEP2)
        print(bold("MOHID Output Comparison  [AUTO MODE]"))
        print(f"  res dir    : {path_a}")
        print(f"  pairs found: {len(pairs)}")
        print(f"  tolerance  : abs={args.tol_abs:.0e}  rel={args.tol_rel:.0e}  (FAIL if BOTH exceeded)")
        print(SEP2)
    else:
        path_b = Path(args.path_b)
        if not path_b.exists():
            print(f"Error: path not found: {path_b}", file=sys.stderr)
            sys.exit(2)
        pairs = [(path_a, path_b, path_b.name)]
        print(SEP2)
        print(bold("MOHID Output Comparison  [EXPLICIT MODE]"))
        print(f"  baseline   : {path_a}")
        print(f"  new        : {path_b}")
        print(f"  tolerance  : abs={args.tol_abs:.0e}  rel={args.tol_rel:.0e}  (FAIL if BOTH exceeded)")
        print(SEP2)

    # ── Run comparisons ───────────────────────────────────────────────────────

    all_results: List[FileResult] = []
    for baseline, new, label in pairs:
        results = _compare_pair(baseline, new, args.tol_abs, args.tol_rel, label)
        for fr in results:
            all_results.append(fr)
            _print_file_result(fr, args.verbose, args.tol_abs, args.tol_rel)

    _print_summary(all_results)

    any_fail = any(r.status in ("fail", "error") for r in all_results)
    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()
