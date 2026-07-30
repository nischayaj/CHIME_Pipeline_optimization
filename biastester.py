"""
check_subset_bias_ks_test.py

Purpose: check whether the 200-file subset used for the pipeline comparison
is a statistically representative sample of the full ~535-file CHIME/FRB
catalog, or whether it differs systematically (selection bias) on DM,
SNR, and/or width.

Method: read ONLY the cheap HDF5 header attributes (no waterfall load, no
detrending/masking/burst search) for every file in the full catalog and
for the subset, then run a two-sample Kolmogorov-Smirnov (KS) test per
field. This does not touch, call, or modify any function in the main
pipeline script -- it is a standalone, read-only diagnostic.

Usage:
    1. Edit DATA_DIR / FILE_PATTERN below if they differ from your setup.
    2. Edit `get_subset_files()` to match EXACTLY how your 200 files were
       chosen (the default assumes "first 200 in sorted filename order" --
       change this if you used a different selection rule).
    3. Run: python check_subset_bias_ks_test.py
"""

import os
import glob
import numpy as np
import h5py
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# CONFIG -- match these to your pipeline's constants
# ---------------------------------------------------------------------------
DATA_DIR = r"C:\Users\nischaya jajodia\Desktop\frbunzipped"  # directory containing the full catalog of HDF5 files
FILE_PATTERN = "*.h5"
N_SUBSET = 200          # how many files your subset used
OUTPUT_DIR = None       # None -> current working directory

# Candidate attribute names to test. Not every CHIME/FRB release uses the
# exact same key names, so we check what's actually present in the files
# before testing anything (see print_available_attrs()).
CANDIDATE_FIELDS = ["dm", "snr", "width_fitb", "width", "fluence"]


# ---------------------------------------------------------------------------
def get_all_files(data_dir=DATA_DIR, pattern=FILE_PATTERN):
    return sorted(glob.glob(os.path.join(data_dir, pattern)))


def get_subset_files(all_files, n_subset=N_SUBSET):
    """[EDIT ME IF NEEDED] Default assumption: subset = first N files in
    sorted filename order (matches "first 200 files"). If your subset was
    chosen a different way (e.g. a specific list, a different sort key,
    or files processed in directory-listing order rather than sorted
    order), change this function to reproduce that exact selection --
    otherwise this diagnostic will be testing the wrong subset."""
    return all_files[:n_subset]


def print_available_attrs(file_name):
    """Prints every attribute actually present on one file, so you can
    confirm the CANDIDATE_FIELDS names above match your data release."""
    with h5py.File(file_name, "r") as f:
        keys = list(f["frb"].attrs.keys())
    print(f"Attributes found in {os.path.basename(file_name)}: {keys}\n")
    return keys


def get_attr_value(file_name, key):
    """Reads a single scalar attribute from one file's header. Returns
    np.nan if the key isn't present or can't be parsed as a scalar --
    never raises, so one malformed file doesn't kill the whole loop."""
    try:
        with h5py.File(file_name, "r") as f:
            attrs = f["frb"].attrs
            if key not in attrs:
                return np.nan
            val = attrs[key][()] if hasattr(attrs[key], "__getitem__") else attrs[key]
            return float(val)
    except Exception:
        return np.nan


def collect_field(files, key):
    return np.array([get_attr_value(fn, key) for fn in files], dtype=float)


def run_ks_comparison(all_files, subset_files, fields):
    """Runs a KS test per field, comparing subset vs full catalog.
    Returns a list of dicts with the results (also printed to console)."""
    results = []
    n_tested = 0
    # first pass: figure out which fields actually have usable data
    usable_fields = []
    for key in fields:
        sample_vals = collect_field(all_files[: min(20, len(all_files))], key)
        if np.isfinite(sample_vals).any():
            usable_fields.append(key)
    if not usable_fields:
        print("No usable fields found among CANDIDATE_FIELDS -- run "
              "print_available_attrs() on a sample file and update "
              "CANDIDATE_FIELDS with the correct key names.")
        return results

    n_tested = len(usable_fields)
    bonferroni_alpha = 0.05 / n_tested

    print(f"Testing {n_tested} field(s): {usable_fields}")
    print(f"Bonferroni-corrected significance threshold: {bonferroni_alpha:.4f} "
          f"(0.05 / {n_tested} tests)\n")

    for key in usable_fields:
        vals_all = collect_field(all_files, key)
        vals_sub = collect_field(subset_files, key)
        vals_all = vals_all[np.isfinite(vals_all)]
        vals_sub = vals_sub[np.isfinite(vals_sub)]

        if len(vals_all) < 5 or len(vals_sub) < 5:
            print(f"[{key}] Skipped -- not enough valid values "
                  f"(all={len(vals_all)}, subset={len(vals_sub)})")
            continue

        stat, pval = ks_2samp(vals_sub, vals_all)
        verdict = "NOT significantly different" if pval > bonferroni_alpha else \
                  "SIGNIFICANTLY DIFFERENT"

        print(f"[{key}] n_all={len(vals_all)}, n_subset={len(vals_sub)}  "
              f"median_all={np.median(vals_all):.3f}  "
              f"median_subset={np.median(vals_sub):.3f}  "
              f"KS_stat={stat:.4f}  p={pval:.4f}  -> {verdict}")

        results.append(dict(
            field=key, n_all=len(vals_all), n_subset=len(vals_sub),
            median_all=float(np.median(vals_all)),
            median_subset=float(np.median(vals_sub)),
            ks_stat=float(stat), p_value=float(pval),
            bonferroni_alpha=bonferroni_alpha,
            significant_after_correction=bool(pval <= bonferroni_alpha),
        ))
    return results


def plot_cdf_comparison(all_files, subset_files, key, save_dir=None):
    """Overlaid empirical CDF plot for one field -- full catalog vs subset.
    Two closely-overlapping curves is the visual you want for the poster."""
    vals_all = collect_field(all_files, key)
    vals_sub = collect_field(subset_files, key)
    vals_all = np.sort(vals_all[np.isfinite(vals_all)])
    vals_sub = np.sort(vals_sub[np.isfinite(vals_sub)])
    if len(vals_all) < 5 or len(vals_sub) < 5:
        return None

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(vals_all, np.arange(1, len(vals_all) + 1) / len(vals_all),
            label=f"Full catalog (n={len(vals_all)})", linewidth=1.6)
    ax.plot(vals_sub, np.arange(1, len(vals_sub) + 1) / len(vals_sub),
            label=f"Subset (n={len(vals_sub)})", linewidth=1.6, linestyle="--")
    ax.set_xlabel(key)
    ax.set_ylabel("Cumulative fraction")
    ax.set_title(f"Empirical CDF: {key}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if save_dir:
        fig.savefig(os.path.join(save_dir, f"ks_cdf_{key}.png"), dpi=150)
    plt.show()
    return fig


if __name__ == "__main__":
    all_files = get_all_files()
    subset_files = get_subset_files(all_files)
    print(f"Full catalog: {len(all_files)} files. Subset: {len(subset_files)} files.\n")

    if not all_files:
        raise SystemExit(f"No files found in {DATA_DIR} matching {FILE_PATTERN} -- "
                          f"check DATA_DIR/FILE_PATTERN before running.")

    # Step 1: confirm the real attribute names on one file
    print_available_attrs(all_files[0])

    # Step 2: run the KS comparison across all usable candidate fields
    results = run_ks_comparison(all_files, subset_files, CANDIDATE_FIELDS)

    # Step 3: save a plot for each field that was actually tested
    for r in results:
        plot_cdf_comparison(all_files, subset_files, r["field"], save_dir=OUTPUT_DIR)

    # Step 4: write results to a small CSV for the appendix/limitations section
    if results:
        import csv
        out_path = os.path.join(OUTPUT_DIR or os.getcwd(), "subset_bias_ks_results.csv")
        with open(out_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to: {out_path}")
