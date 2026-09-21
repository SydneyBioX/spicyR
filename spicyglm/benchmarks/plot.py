"""Plot benchmark results from benchmarks/results/timings.csv.

Usage: python benchmarks/plot.py benchmarks/results/timings.csv benchmarks/results

Writes benchmarks.pdf / .png (the two-panel headline figure), benchmarks_detailed.pdf
(the four detailed figures as panels A-D), and one PNG per detailed figure. Points are medians over repeats; whiskers span the
fastest and slowest repeat. Speed-ups compare like with like: 1 R core against
1 spicyglm thread, and 4 R cores against 4 threads.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
# colour = implementation (categorical slots 1-2, validated all-pairs); line style + marker = workers
COLOR = {"R spicyGLM()": "#2a78d6", "spicyglm": "#eb6834"}
STYLE = {1: dict(linestyle="-", marker="o"), 4: dict(linestyle="--", marker="^")}
WORKERS = {("R spicyGLM()", 1): "1 core", ("R spicyGLM()", 4): "4 cores", ("spicyglm", 1): "1 thread",
           ("spicyglm", 4): "4 threads"}
SERIES = [("R spicyGLM()", 1), ("R spicyGLM()", 4), ("spicyglm", 1), ("spicyglm", 4)]
MEM_LIMIT_GB = 6  # must match benchmark.py
FAMILIES = ["poisson", "binomial"]
FAMILY_TITLE = {"poisson": "Poisson (fixed radius)", "binomial": "Binomial (k nearest neighbours)"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK_2, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
    "axes.titleweight": "bold", "legend.frameon": False, "legend.handlelength": 3.2,
    "pdf.fonttype": 42,  # editable text in the PDF
})
CELLS = FuncFormatter(lambda v, _: f"{v / 1e6:g}M" if v >= 1e6 else f"{v / 1e3:g}k")


def secs(t):
    return f"{t:.2g} s" if t < 10 else f"{t:.0f} s"


def label(impl, workers):
    return f"{impl}, {WORKERS[(impl, workers)]}"


def series(sub, impl, workers):
    return sub[(sub.implementation == impl) & (sub.workers == workers) & ~sub.stopped_for_memory].sort_values("n_cells")


def end_label(ax, x, y, text):
    ax.annotate(text, (np.asarray(x)[-1], np.asarray(y)[-1]), xytext=(8, 0), textcoords="offset points",
                va="center", fontsize=9, color=INK_2)


def stopped_note(df):
    if not df.stopped_for_memory.any():
        return None
    ok = df[(df.implementation == "R spicyGLM()") & (df.workers == 4) & ~df.stopped_for_memory]
    return (f"R spicyGLM() with 4 cores was stopped at the {MEM_LIMIT_GB:g} GB memory limit above "
            f"{ok.n_cells.max() / 1e3:g}k cells and on Schürch 2020 (memory summed over its forked workers).")


def cell_axes(axes, ylabel):
    for ax, family in zip(axes, FAMILIES):
        ax.set_title(FAMILY_TITLE[family])
        ax.set_xlabel("Cells (all 36 cell-type pairs)")
        ax.xaxis.set_major_formatter(CELLS)
    axes[0].set_ylabel(ylabel)
    axes[0].legend(loc="upper left")


# ---- panels: each draws into a Figure or SubFigure -------------------------

def draw_timings(fig, synth, title):
    axes = fig.subplots(1, 2, sharey=True)
    for ax, family in zip(axes, FAMILIES):
        sub = synth[synth.family == family]
        for impl, workers in SERIES:
            s = series(sub, impl, workers)
            ax.plot(s.n_cells, s.seconds_median, color=COLOR[impl], linewidth=2, markersize=7,
                    markeredgecolor=SURFACE, markeredgewidth=1.5, label=label(impl, workers), **STYLE[workers])
            ax.vlines(s.n_cells, s.seconds_min, s.seconds_max, color=COLOR[impl], linewidth=1.2, alpha=0.8)
            end_label(ax, s.n_cells, s.seconds_median, secs(s.seconds_median.iloc[-1]))
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(sub.n_cells.min() / 1.4, sub.n_cells.max() * 2.4)
    cell_axes(axes, "Fitting time (s, median, log scale)")
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold")


def draw_speedup(fig, synth, title):
    axes = fig.subplots(1, 2, sharey=True)
    top = 0
    for ax, family in zip(axes, FAMILIES):
        sub = synth[synth.family == family]
        for workers in (1, 4):
            r = series(sub, "R spicyGLM()", workers).set_index("n_cells").seconds_median
            p = series(sub, "spicyglm", workers).set_index("n_cells").seconds_median
            common = r.index.intersection(p.index)
            speed = (r.loc[common] / p.loc[common]).to_numpy()
            top = max(top, speed.max())
            ax.plot(common, speed, color=COLOR["spicyglm"], linewidth=2, markersize=7, markeredgecolor=SURFACE,
                    markeredgewidth=1.5, label=f"{WORKERS[('R spicyGLM()', workers)]} vs {WORKERS[('spicyglm', workers)]}",
                    **STYLE[workers])
            end_label(ax, common, speed, f"{speed[-1]:.0f}x")
        ax.set_xscale("log")
        ax.set_xlim(sub.n_cells.min() / 1.4, sub.n_cells.max() * 2.2)
    axes[0].set_ylim(0, top * 1.12)
    cell_axes(axes, "Speed-up (R time / spicyglm time)")
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold")


def draw_memory(fig, synth, title):
    axes = fig.subplots(1, 2, sharey=True)
    for ax, family in zip(axes, FAMILIES):
        sub = synth[synth.family == family]
        for impl, workers in SERIES:
            s = series(sub, impl, workers)
            gb = s.peak_rss_mb / 1024
            ax.plot(s.n_cells, gb, color=COLOR[impl], linewidth=2, markersize=7, markeredgecolor=SURFACE,
                    markeredgewidth=1.5, label=label(impl, workers), **STYLE[workers])
            if impl.startswith("R") or workers == 4:  # the two spicyglm lines overlap: label once
                end_label(ax, s.n_cells, gb, f"{gb.iloc[-1]:.2g} GB")
        ax.set_xscale("log")
        ax.set_xlim(sub.n_cells.min() / 1.4, sub.n_cells.max() * 2.4)
    axes[0].set_ylim(0, synth.peak_rss_mb.max() / 1024 * 1.12)
    cell_axes(axes, "Peak memory (GB)")
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold")


def draw_real(fig, real, title):
    axes = fig.subplots(1, 2, sharex=True)
    order = SERIES[::-1]
    for ax, family in zip(axes, FAMILIES):
        sub = real[real.family == family]
        for i, (impl, workers) in enumerate(order):
            row = sub[(sub.implementation == impl) & (sub.workers == workers)].iloc[0]
            if row.stopped_for_memory:
                ax.annotate(f"stopped at the {MEM_LIMIT_GB:g} GB memory limit", (real.seconds_min.min(), i),
                            xytext=(4, 0), textcoords="offset points", va="center", fontsize=9, color=INK_2)
                continue
            r_same = sub[(sub.implementation == "R spicyGLM()") & (sub.workers == workers)].iloc[0]
            ax.hlines(i, row.seconds_min, row.seconds_max, color=COLOR[impl], linewidth=1.2)
            ax.plot(row.seconds_median, i, color=COLOR[impl], markersize=10, markeredgecolor=SURFACE,
                    markeredgewidth=1.5, linestyle="none", marker=STYLE[workers]["marker"], zorder=3)
            text = secs(row.seconds_median)
            if impl == "spicyglm" and not r_same.stopped_for_memory:
                text += f"  ({r_same.seconds_median / row.seconds_median:.0f}x vs R, {WORKERS[('R spicyGLM()', workers)]})"
            ax.annotate(text, (row.seconds_max, i), xytext=(10, 0), textcoords="offset points", va="center",
                        fontsize=9, color=INK_2)
        ax.set_yticks(range(len(order)), [label(*s) for s in order], color=INK_2)
        ax.set(xscale="log", xlabel="Fitting time (s, median, log scale)", ylim=(-0.6, len(order) - 0.4))
        ax.set_title(FAMILY_TITLE[family] + (", with diagnostics" if sub.diagnostics.iloc[0] else ""))
        ax.grid(axis="y", visible=False)
    axes[0].set_xlim(real.seconds_min.min() / 3, real.seconds_max.max() * 60)
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold")


def draw_summary(fig, df):
    rows = [("synthetic_250000", "poisson", "Synthetic, 250k cells, Poisson"),
            ("synthetic_250000", "binomial", "Synthetic, 250k cells, Binomial"),
            ("synthetic_1000000", "poisson", "Synthetic, 1M cells, Poisson"),
            ("synthetic_1000000", "binomial", "Synthetic, 1M cells, Binomial"),
            ("Schurch 2020", "poisson", "Schürch 2020 (real), Poisson + diagnostics"),
            ("Schurch 2020", "binomial", "Schürch 2020 (real), Binomial")]
    one = df[df.workers == 1].set_index(["dataset", "family", "implementation"]).seconds_median
    ax = fig.subplots()
    for i, (dataset, family, name) in enumerate(rows[::-1]):
        r, p = one[(dataset, family, "R spicyGLM()")], one[(dataset, family, "spicyglm")]
        ax.barh(i, r / p, height=0.62, color=COLOR["spicyglm"])
        ax.annotate(f"{r / p:.0f}x   (R {secs(r)} → {secs(p)})", (r / p, i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9.5, color=INK_2)
    ax.set_yticks(range(len(rows)), [name for *_, name in rows[::-1]], color=INK_2)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, ax.get_xlim()[1] * 1.45)
    ax.set_xlabel("Speed-up over R spicyGLM() (median fitting time, 1 core each)")
    fig.suptitle("spicyglm (C++ core) vs R spicyGLM()", x=0.01, ha="left", fontweight="bold")


# ---- output ------------------------------------------------------------------

def draw_story(fig, df):
    """Two-panel headline figure: 1 core each, scaling on synthetic data and the real dataset."""
    one = df[(df.workers == 1) & ~df.stopped_for_memory]
    synth = one[one.dataset.str.startswith("synthetic")]
    real = one[~one.dataset.str.startswith("synthetic")]
    fam_style = {"poisson": dict(linestyle="-", marker="o"), "binomial": dict(linestyle="--", marker="s")}
    fam_name = {"poisson": "Poisson", "binomial": "Binomial"}
    impl_name = {"R spicyGLM()": "R spicyGLM()", "spicyglm": "spicyglm"}

    ax_a, ax_b = fig.subplots(1, 2, width_ratios=[1.25, 1])

    # A: fitting time against dataset size
    ends = {}
    for impl in ("R spicyGLM()", "spicyglm"):
        for family in FAMILIES:
            s = synth[(synth.implementation == impl) & (synth.family == family)].sort_values("n_cells")
            ax_a.plot(s.n_cells, s.seconds_median, color=COLOR[impl], linewidth=2.2, markersize=6.5,
                      markeredgecolor=SURFACE, markeredgewidth=1.3, **fam_style[family])
            ends[(impl, family)] = (s.n_cells.iloc[-1], s.seconds_median.iloc[-1])
    for (impl, family), (x, y) in ends.items():
        # nudge the two labels of an implementation apart so they never collide
        dy = {"poisson": -7, "binomial": 7} if impl == "spicyglm" else {"poisson": -6, "binomial": 6}
        ax_a.annotate(f"{impl_name[impl]}, {fam_name[family]}  {secs(y)}", (x, y), xytext=(9, dy[family]),
                      textcoords="offset points", va="center", fontsize=9, color=COLOR[impl], fontweight="bold")
    x_end = synth.n_cells.max()
    lo = min(ends[("spicyglm", f)][1] for f in FAMILIES)
    hi = max(ends[("R spicyGLM()", f)][1] for f in FAMILIES)
    speed = [ends[("R spicyGLM()", f)][1] / ends[("spicyglm", f)][1] for f in FAMILIES]
    ax_a.annotate("", xy=(x_end * 0.86, lo * 1.25), xytext=(x_end * 0.86, hi / 1.25),
                  arrowprops=dict(arrowstyle="<->", color=INK_2, linewidth=1.2))
    ax_a.annotate(f"{min(speed):.0f}-{max(speed):.0f}x\nfaster", (x_end * 0.8, np.sqrt(lo * hi)), ha="right",
                  va="center", fontsize=11, fontweight="bold", color=INK)
    ax_a.set(xscale="log", yscale="log", xlabel="Cells in the dataset (all 36 cell-type pairs fitted)",
             ylabel="Fitting time (seconds, log scale)")
    ax_a.set_xlim(synth.n_cells.min() / 1.3, x_end * 5.5)
    ax_a.xaxis.set_major_formatter(CELLS)
    ax_a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax_a.set_title("A   Faster at every dataset size (synthetic data)")

    # B: real dataset, R and spicyglm joined per model
    rows = [("poisson", "Poisson\n+ diagnostics" if real[real.family == "poisson"].diagnostics.iloc[0] else "Poisson"),
            ("binomial", "Binomial")]
    for i, (family, name) in enumerate(rows[::-1]):
        sub = real[real.family == family].set_index("implementation").seconds_median
        r, p = sub["R spicyGLM()"], sub["spicyglm"]
        ax_b.plot([p, r], [i, i], color=AXIS, linewidth=2.5, zorder=1)
        ax_b.plot(r, i, "o", color=COLOR["R spicyGLM()"], markersize=11, markeredgecolor=SURFACE, zorder=3)
        ax_b.plot(p, i, "o", color=COLOR["spicyglm"], markersize=11, markeredgecolor=SURFACE, zorder=3)
        ax_b.annotate(secs(r), (r, i), xytext=(0, -15), textcoords="offset points", ha="center", va="top",
                      fontsize=9, color=COLOR["R spicyGLM()"], fontweight="bold")
        ax_b.annotate(secs(p), (p, i), xytext=(0, -15), textcoords="offset points", ha="center", va="top",
                      fontsize=9, color=COLOR["spicyglm"], fontweight="bold")
        ax_b.annotate(f"{r / p:.0f}x faster", (np.sqrt(r * p), i), xytext=(0, 9), textcoords="offset points",
                      ha="center", va="bottom", fontsize=11, fontweight="bold", color=INK)
    ax_b.set_yticks(range(len(rows)), [name for _, name in rows[::-1]], color=INK_2, fontsize=10)
    ax_b.set(xscale="log", xlabel="Fitting time (seconds, log scale)", ylim=(-0.7, len(rows) - 0.4))
    ax_b.set_xlim(real.seconds_median.min() / 3, real.seconds_median.max() * 3)
    ax_b.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax_b.grid(axis="y", visible=False)
    n_real = f"{df[~df.dataset.str.startswith('synthetic')].n_cells.iloc[0] / 1e3:.0f}k"
    ax_b.set_title(f"B   Real data: Schürch 2020 CRC ({n_real} cells)")

    real_speed = [real[(real.family == f) & (real.implementation == "R spicyGLM()")].seconds_median.iloc[0]
                  / real[(real.family == f) & (real.implementation == "spicyglm")].seconds_median.iloc[0]
                  for f in FAMILIES]
    lo_all, hi_all = min(real_speed + speed), max(real_speed + speed)
    fig.suptitle(f"spicyglm fits the same model {lo_all:.0f}-{hi_all:.0f}x faster than R spicyGLM()",
                 fontsize=14, fontweight="bold")
    big = df[(df.n_cells == df.n_cells.max()) & (df.workers == 1)].groupby("implementation").peak_rss_mb.max() / 1024
    fig.supxlabel(
        f"Both on 1 core; median fitting time over repeated runs (5 spicyglm, 3 R), estimates match R. "
        f"At {df.n_cells.max() / 1e6:g}M cells spicyglm also used {big['spicyglm']:.1f} GB of memory vs "
        f"{big['R spicyGLM()']:.1f} GB for R.\nOne laptop (Apple silicon, 10 cores, 24 GB).",
        fontsize=8.5, color=INK_2, multialignment="center")


def main():
    df = pd.read_csv(sys.argv[1])
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    synth = df[df.dataset.str.startswith("synthetic")]
    real = df[~df.dataset.str.startswith("synthetic")]
    note = stopped_note(df)
    n_real = f"{real.n_cells.iloc[0]:,}" if not real.empty else ""

    panels = [
        ("timings", 4.6, lambda f, t: draw_timings(f, synth, t),
         "Fitting time: R spicyGLM() vs spicyglm (synthetic data, 5,000 cells per image)"),
        ("speedup", 4.4, lambda f, t: draw_speedup(f, synth, t),
         "Speed-up of spicyglm over R spicyGLM(), same number of workers"),
        ("memory", 4.4, lambda f, t: draw_memory(f, synth, t),
         "Peak memory: whole process tree, including start-up and reading the CSV"),
    ]
    if not real.empty:
        panels.append(("schurch", 3.9, lambda f, t: draw_real(f, real, t),
                       f"Schürch 2020 colorectal cancer CODEX data ({n_real} cells, 140 images, 435 pairs)"))

    for name, height, draw, title in panels:
        fig = plt.figure(figsize=(11.5, height + (0.35 if note else 0)), layout="constrained")
        draw(fig, title)
        if note:
            fig.supxlabel(note, x=0.01, ha="left", fontsize=8.5, color=INK_2)
        fig.savefig(out / f"{name}.png", dpi=160)
        plt.close(fig)

    fig = plt.figure(figsize=(9.5, 4.2), layout="constrained")
    draw_summary(fig, df)
    fig.savefig(out / "summary.png", dpi=160)
    plt.close(fig)

    # all detailed panels on one page
    heights = [h for _, h, _, _ in panels]
    fig = plt.figure(figsize=(11.5, sum(heights) + 0.9), layout="constrained")
    fig.get_layout_engine().set(hspace=0.04)
    subfigs = fig.subfigures(len(panels), 1, height_ratios=heights)
    for letter, sub, (_, _, draw, title) in zip("ABCD", subfigs, panels):
        draw(sub, f"{letter}   {title}")
    fig.supxlabel("Medians over repeated runs (5 spicyglm, 3 R); whiskers show the fastest and slowest run. "
                  "Fitting time only, one laptop (Apple silicon, 10 cores, 24 GB).\n" + (note or ""),
                  x=0.01, ha="left", fontsize=8.5, color=INK_2)
    fig.savefig(out / "benchmarks_detailed.pdf")
    plt.close(fig)

    # headline figure
    fig = plt.figure(figsize=(12, 5.2), layout="constrained")
    draw_story(fig, df)
    fig.savefig(out / "benchmarks.pdf")
    fig.savefig(out / "benchmarks.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
