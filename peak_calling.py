"""Toy CUT&Tag analysis: background-threshold peak calling and depth normalization.

Part 1 simulates a target (antibody) sample and an IgG control, sets a threshold
from the control's background distribution (a simplified version of the idea
behind SEACR), and measures how well true binding sites are recovered.

Part 2 shows why sequencing depth must be normalized (CPM) before comparing
samples, and computes FRiP (fraction of reads in peaks) as a quality metric.
"""
import numpy as np
import matplotlib.pyplot as plt


def call_peaks(seed=42, n_bins=2000, n_peaks=50):
    rng = np.random.default_rng(seed)
    true_peaks = rng.choice(n_bins, size=n_peaks, replace=False)
    is_peak = np.zeros(n_bins, dtype=bool)
    is_peak[true_peaks] = True

    # Target sample: low background everywhere, high counts at real binding sites
    target = rng.poisson(2, n_bins)
    target[is_peak] = rng.poisson(15, is_peak.sum())
    # IgG control: background only
    igg = rng.poisson(2, n_bins)

    # Threshold from the control's background distribution
    threshold = np.percentile(igg, 99.9)
    called = target > threshold

    tp = (called & is_peak).sum()
    fp = (called & ~is_peak).sum()
    fn = (~called & is_peak).sum()
    print("Part 1: peak calling")
    print(f"  threshold = {threshold:.1f} reads")
    print(f"  called peaks = {called.sum()}, true positives = {tp}, "
          f"false positives = {fp}, missed = {fn}")
    print(f"  precision = {tp / (tp + fp):.2f}, recall = {tp / (tp + fn):.2f}")

    plt.figure(figsize=(9, 3))
    plt.plot(target, lw=0.6, color="gray", label="target sample")
    plt.scatter(np.where(called)[0], target[called], s=12, color="red", label="called peaks")
    plt.axhline(threshold, ls="--", color="blue", label=f"threshold = {threshold:.0f}")
    plt.xlabel("genomic bin")
    plt.ylabel("reads")
    plt.legend(loc="upper right")
    plt.title("Background-threshold peak calling")
    plt.savefig("peaks.png", dpi=150, bbox_inches="tight")
    print("  saved peaks.png")


def normalize(seed=0, n_bins=2000, n_peaks=50):
    rng = np.random.default_rng(seed)
    peaks = np.zeros(n_bins, dtype=bool)
    peaks[rng.choice(n_bins, n_peaks, replace=False)] = True

    def sample(depth_scale):
        c = rng.poisson(2 * depth_scale, n_bins)
        c[peaks] = rng.poisson(15 * depth_scale, peaks.sum())
        return c

    def cpm(c):
        return c / c.sum() * 1e6

    def frip(c):
        return c[peaks].sum() / c.sum()

    a = sample(1.0)  # sample A
    b = sample(3.0)  # sample B: same biology, sequenced 3x deeper
    print("Part 2: normalization")
    print(f"  raw totals: A = {a.sum()}, B = {b.sum()}")
    print(f"  mean raw count in peaks: A = {a[peaks].mean():.1f}, B = {b[peaks].mean():.1f}")
    print(f"  mean CPM in peaks:       A = {cpm(a)[peaks].mean():.0f}, B = {cpm(b)[peaks].mean():.0f}")
    print(f"  FRiP: A = {frip(a):.2f}, B = {frip(b):.2f}")


if __name__ == "__main__":
    call_peaks()
    normalize()
