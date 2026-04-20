#!/usr/bin/env python3
"""
scrape_versions.py
------------------
Extract exact installed versions of every tool used by the gcl_illumina_qc
Nextflow pipeline by:

  1. Parsing nextflow.config for every `conda = '...'` directive and the
     label it belongs to.
  2. Walking the Nextflow conda cache directory and reading the
     conda-meta/*.json metadata in each environment to get actual
     installed versions (not just spec versions, which may be loose
     like "bioconda::samtools" with no version pin).
  3. Matching envs to labels by package overlap (Nextflow names env
     directories by a hash of the spec + channels, so we can't go by
     name directly).
  4. Printing (a) a per-label version table and (b) a citation-ready
     summary of the tools you actually reference in the methods section.

Usage
-----
    python scrape_versions.py \\
        --config /path/to/nextflow.config \\
        --cache  /work/birdlab/.conda_builds

If you leave --cache off, it will be read from the `conda.cacheDir`
line in nextflow.config.

Output goes to stdout; redirect to a file if you want to save it:

    python scrape_versions.py --config nextflow.config > versions.md
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# 1. Parse nextflow.config
# ---------------------------------------------------------------------------

# Matches a `withLabel: name { ... }` block loosely enough to find its conda
# line without needing a full Groovy parser. We assume the closing brace for
# a label block sits on its own line (true in the provided config).
LABEL_RE = re.compile(r"withLabel:\s*['\"]?([A-Za-z0-9_]+)['\"]?\s*\{")
CONDA_RE = re.compile(r"conda\s*=\s*['\"]([^'\"]+)['\"]")
CACHE_RE = re.compile(r"cacheDir\s*=\s*['\"]([^'\"]+)['\"]")


def parse_config(config_path: Path) -> Tuple[Dict[str, List[str]], Optional[str]]:
    """Return (label -> list of package specs, conda cache dir)."""
    text = config_path.read_text()

    cache_match = CACHE_RE.search(text)
    cache_dir = cache_match.group(1) if cache_match else None

    # Scan line by line, tracking the current label.
    label_specs: Dict[str, List[str]] = OrderedDict()
    current_label: Optional[str] = None
    brace_depth = 0  # depth inside the current label block

    for raw_line in text.splitlines():
        line = raw_line.strip()

        # strip // comments
        if "//" in line:
            line = line.split("//", 1)[0].strip()
        if not line:
            continue

        label_hit = LABEL_RE.search(line)
        if label_hit and current_label is None:
            current_label = label_hit.group(1)
            brace_depth = line.count("{") - line.count("}")
            continue

        if current_label is not None:
            brace_depth += line.count("{")
            brace_depth -= line.count("}")

            conda_hit = CONDA_RE.search(line)
            if conda_hit:
                # split on whitespace; each token is like "bioconda::fastp" or
                # "conda-forge::r-base=4.3"
                specs = conda_hit.group(1).split()
                label_specs[current_label] = specs

            if brace_depth <= 0:
                current_label = None
                brace_depth = 0

    return label_specs, cache_dir


def pkg_name_from_spec(spec: str) -> str:
    """'bioconda::r-base=4.3' -> 'r-base'."""
    # drop channel prefix
    if "::" in spec:
        spec = spec.split("::", 1)[1]
    # drop version constraint
    for sep in ("=", ">", "<", "!", " "):
        if sep in spec:
            spec = spec.split(sep, 1)[0]
    return spec.strip().lower()


# ---------------------------------------------------------------------------
# 2. Walk the conda cache and read installed versions
# ---------------------------------------------------------------------------

def read_env_packages(env_dir: Path) -> Dict[str, Dict[str, str]]:
    """Return {pkg_name: {version, build, channel}} for one env."""
    meta = env_dir / "conda-meta"
    if not meta.is_dir():
        return {}

    pkgs: Dict[str, Dict[str, str]] = {}
    for j in meta.glob("*.json"):
        try:
            data = json.loads(j.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        name = data.get("name")
        if not name:
            continue
        pkgs[name.lower()] = {
            "version": data.get("version", "?"),
            "build": data.get("build", ""),
            "channel": data.get("channel", ""),
        }
    return pkgs


def scan_cache(cache_dir: Path) -> Dict[Path, Dict[str, Dict[str, str]]]:
    """Return {env_dir: {pkg: {version,...}}} for every env in cache."""
    envs: Dict[Path, Dict[str, Dict[str, str]]] = {}
    if not cache_dir.is_dir():
        return envs
    for child in sorted(cache_dir.iterdir()):
        if not child.is_dir():
            continue
        pkgs = read_env_packages(child)
        if pkgs:
            envs[child] = pkgs
    return envs


# ---------------------------------------------------------------------------
# 3. Match label -> env by package overlap (best Jaccard on primary pkgs)
# ---------------------------------------------------------------------------

def best_env_for_label(
    specs: List[str],
    envs: Dict[Path, Dict[str, Dict[str, str]]],
) -> Optional[Path]:
    wanted = {pkg_name_from_spec(s) for s in specs}
    if not wanted:
        return None

    best: Tuple[float, Optional[Path]] = (0.0, None)
    for env_dir, pkgs in envs.items():
        installed = set(pkgs.keys())
        # require all wanted packages to be present in the env
        if not wanted.issubset(installed):
            continue
        # among envs that satisfy the label, prefer the smallest
        # (most specific) environment to avoid matching a fat shared env
        score = len(wanted) / len(installed)
        if score > best[0]:
            best = (score, env_dir)
    return best[1]


# ---------------------------------------------------------------------------
# 4. Report
# ---------------------------------------------------------------------------

# Tools that appear in the methods section — script will print a tidy
# citation-ready block just for these.
CITATION_TOOLS = [
    # preprocessing
    "fastp",
    "bbmap",          # BBTools (clumpify, repair)
    "fastq-screen",
    "bowtie2",
    # alignment / BAM
    "bwa-mem2",
    "samtools",
    # QC aggregation
    "fastqc",
    "multiqc",
    # R + key analysis packages for read-retention modelling
    "r-base",
    "r-tidyverse",
    "r-glmmtmb",
    "r-emmeans",
    "r-ggplot2",
    # reporting
    "pandoc",
    "ncbi-datasets-cli",
]


def fmt_table(rows: List[Tuple[str, ...]], headers: Tuple[str, ...]) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    header_line = "| " + " | ".join(h.ljust(w) for h, w in zip(headers, widths)) + " |"
    body = [
        "| " + " | ".join(str(c).ljust(w) for c, w in zip(row, widths)) + " |"
        for row in rows
    ]
    return "\n".join([header_line, sep, *body])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=Path("nextflow.config"),
                    help="Path to nextflow.config (default: ./nextflow.config)")
    ap.add_argument("--cache", type=Path, default=None,
                    help="Conda cache dir. Default: read from nextflow.config.")
    args = ap.parse_args()

    if not args.config.is_file():
        print(f"ERROR: config not found: {args.config}", file=sys.stderr)
        return 1

    label_specs, cfg_cache = parse_config(args.config)
    cache_dir = args.cache or (Path(cfg_cache) if cfg_cache else None)

    print(f"# Versions for {args.config}\n")
    print(f"- Config: `{args.config}`")
    print(f"- Conda cache: `{cache_dir}`")
    print(f"- Labels found: {len(label_specs)}\n")

    if cache_dir is None or not cache_dir.is_dir():
        print(f"ERROR: conda cache dir not found or not readable: {cache_dir}",
              file=sys.stderr)
        print("\nLabels and their declared specs (no envs to resolve against):\n")
        for label, specs in label_specs.items():
            print(f"- **{label}**: {' '.join(specs)}")
        return 2

    envs = scan_cache(cache_dir)
    print(f"- Envs in cache: {len(envs)}\n")

    # ---- Per-label resolved versions -------------------------------------
    print("## Per-label resolved versions\n")

    # Aggregate pkg -> set of (version, build) across matched envs
    tool_versions: Dict[str, Dict[str, set]] = {}

    for label, specs in label_specs.items():
        env_dir = best_env_for_label(specs, envs)
        print(f"### `{label}`\n")
        print(f"Declared: `{' '.join(specs)}`\n")
        if env_dir is None:
            print("_No matching conda env found in cache._\n")
            continue
        print(f"Matched env: `{env_dir.name}`\n")

        pkgs = envs[env_dir]
        rows = []
        for spec in specs:
            name = pkg_name_from_spec(spec)
            info = pkgs.get(name)
            if info is None:
                rows.append((name, "MISSING", "", ""))
                continue
            rows.append((name, info["version"], info["build"], info["channel"]))
            tool_versions.setdefault(name, set()).add(info["version"])
        print(fmt_table(rows, ("package", "version", "build", "channel")))
        print()

    # ---- Citation-ready summary ------------------------------------------
    print("\n## Citation-ready summary\n")
    print("Versions of tools referenced in the Methods section:\n")
    rows = []
    for tool in CITATION_TOOLS:
        versions = sorted(tool_versions.get(tool, set()))
        if not versions:
            rows.append((tool, "not found in any matched env"))
        else:
            rows.append((tool, ", ".join(versions)))
    print(fmt_table(rows, ("tool", "version(s)")))
    print()

    # ---- Non-conda dependencies (containers, etc.) -----------------------
    # Flag any labels that use containers (visible to user but not scraped).
    container_labels = []
    for raw_line in args.config.read_text().splitlines():
        m = re.search(r"container\s*=\s*['\"]([^'\"]+)['\"]", raw_line)
        if m:
            container_labels.append(m.group(1))
    if container_labels:
        print("\n## Container-based processes (not scraped)\n")
        for c in container_labels:
            print(f"- `{c}`")

    return 0


if __name__ == "__main__":
    sys.exit(main())