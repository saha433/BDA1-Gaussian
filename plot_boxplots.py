import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # no display needed — saves to PNG
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

C_OVERALL = "#2563eb"   # blue
C_10MIN   = "#16a34a"   # green
C_1MIN    = "#dc2626"   # red
C_GAUSS   = "#7c3aed"   # violet

GAUSS_MU    = 0.5
GAUSS_SIGMA = 0.15

def row_to_boxdict(row):

    q1, med, q3 = float(row["p25"]), float(row["median"]), float(row["p75"])
    iqr = q3 - q1
    lo  = q1 - 1.5 * iqr
    hi  = q3 + 1.5 * iqr
    mn  = float(row["min"])
    mx  = float(row["max"])
    return dict(med=med, q1=q1, q3=q3,
                whislo=max(lo, mn), whishi=min(hi, mx),
                fliers=[], mean=float(row["mean"]))

def load_csv(fname):
    if not os.path.exists(fname):
        print(f"  [WARNING] {fname} not found — run gaussian_stream_stats.exe first.")
        return []
    with open(fname, newline="") as f:
        return list(csv.DictReader(f))

global_rows  = load_csv("global_stats.csv")
minute_rows  = load_csv("minute_stats.csv")
ten_min_rows = load_csv("ten_min_stats.csv")

overall_row   = global_rows[0]  if global_rows  else None
block_rows    = global_rows[1:] if global_rows  else []   

fig = plt.figure(figsize=(22, 18), facecolor="#0f172a")
gs  = GridSpec(3, 2, figure=fig,
               left=0.06, right=0.97, top=0.93, bottom=0.06,
               wspace=0.35, hspace=0.55)

LABEL_KW   = dict(color="#e2e8f0", fontsize=9)
TITLE_KW   = dict(color="#f1f5f9", fontsize=11, fontweight="bold", pad=8)
SPINE_COL  = "#334155"
TICK_COL   = "#94a3b8"
GRID_KW    = dict(color="#1e293b", linestyle="--", linewidth=0.7, alpha=0.7)

def style_ax(ax, title):
    ax.set_facecolor("#1e293b")
    for sp in ax.spines.values():
        sp.set_color(SPINE_COL)
    ax.tick_params(colors=TICK_COL, labelsize=8)
    ax.set_title(title, **TITLE_KW)
    ax.grid(True, axis="y", **GRID_KW)
    ax.set_ylim(0.0, 1.0)

def make_bp(ax, stats_list, color, labels):
    bp = ax.bxp(stats_list, showfliers=True, patch_artist=True,
                medianprops=dict(color="#fbbf24", linewidth=1.8),
                boxprops=dict(facecolor=color+"33", edgecolor=color, linewidth=1.2),
                whiskerprops=dict(color=color, linewidth=1.0, linestyle="--"),
                capprops=dict(color=color, linewidth=1.2),
                flierprops=dict(marker=".", color=color, markersize=2, alpha=0.4),
                showmeans=True,
                meanprops=dict(marker="D", markerfacecolor="#f472b6",
                               markeredgecolor="#f472b6", markersize=5))
    ax.set_xticks(range(1, len(labels)+1))
    ax.set_xticklabels(labels, rotation=45, ha="right",
                       color=TICK_COL, fontsize=7.5)
    ax.axhline(GAUSS_MU, color="#fbbf24", linestyle=":", linewidth=0.9, alpha=0.8)
    return bp

# ── Panel 1: Overall (60-min) box plot ───────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
style_ax(ax1, "Panel A — Overall 60-minute Box Plot")
if overall_row:
    make_bp(ax1, [row_to_boxdict(overall_row)], C_OVERALL, ["60 min\n(Overall)"])
ax1.set_ylabel("Value  (Gaussian, μ=0.5, σ=0.15)", **LABEL_KW)

# ── Panel 2: 10-minute blocks ────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
style_ax(ax2, "Panel B — 10-minute Block Box Plots")
if block_rows:
    bd    = [row_to_boxdict(r) for r in block_rows]
    labs  = [f"Block {i+1}\n(min {i*10+1}-{(i+1)*10})" for i in range(len(block_rows))]
    make_bp(ax2, bd, C_10MIN, labs)
ax2.set_ylabel("Value", **LABEL_KW)

# ── Panel 3: Per-minute (all 60) ─────────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, :])
style_ax(ax3, "Panel C — Per-minute Box Plots (all 60 minutes)")
ax3.set_facecolor("#1e293b")
if minute_rows:
    bd   = [row_to_boxdict(r) for r in minute_rows]
    labs = [f"m{r['minute']}" for r in minute_rows]
    make_bp(ax3, bd, C_1MIN, labs)
ax3.set_xlabel("Minute", color=TICK_COL, fontsize=8)
ax3.set_ylabel("Value", **LABEL_KW)

# ── Panel 4: Gaussian bell curve overlay ─────────────────────────────────────
ax4 = fig.add_subplot(gs[2, 0])
ax4.set_facecolor("#1e293b")
style_ax(ax4, "Panel D — Theoretical Gaussian PDF  (μ=0.5, σ=0.15)")
ax4.set_ylim(None, None)   # auto for PDF

x = np.linspace(0.0, 1.0, 500)
pdf = (1 / (GAUSS_SIGMA * math.sqrt(2*math.pi))) * \
      np.exp(-0.5 * ((x - GAUSS_MU)/GAUSS_SIGMA)**2)

ax4.plot(x, pdf, color=C_GAUSS, linewidth=2.2, label="Gaussian PDF")
ax4.fill_between(x, pdf, alpha=0.15, color=C_GAUSS)

# shade IQR region  (mu ± 0.6745*sigma)
q1_th = GAUSS_MU - 0.6745*GAUSS_SIGMA
q3_th = GAUSS_MU + 0.6745*GAUSS_SIGMA
x_iqr = np.linspace(q1_th, q3_th, 200)
p_iqr = (1/(GAUSS_SIGMA*math.sqrt(2*math.pi))) * \
        np.exp(-0.5*((x_iqr-GAUSS_MU)/GAUSS_SIGMA)**2)
ax4.fill_between(x_iqr, p_iqr, alpha=0.35, color=C_10MIN, label="IQR (50% of data)")

ax4.axvline(GAUSS_MU, color="#fbbf24", linestyle="--", linewidth=1.2, label="μ = 0.5")
ax4.axvline(q1_th,    color=C_1MIN,   linestyle=":",  linewidth=1.0, label="Q1 / Q3")
ax4.axvline(q3_th,    color=C_1MIN,   linestyle=":",  linewidth=1.0)

iqr_th = q3_th - q1_th
lo_th  = q1_th - 1.5*iqr_th
hi_th  = q3_th + 1.5*iqr_th
ax4.axvline(max(lo_th,0), color="#f97316", linestyle="-.", linewidth=1.0, label="Tukey fences")
ax4.axvline(min(hi_th,1), color="#f97316", linestyle="-.", linewidth=1.0)

ax4.set_xlabel("Value", color=TICK_COL, fontsize=8)
ax4.set_ylabel("Probability Density", **LABEL_KW)
leg = ax4.legend(fontsize=7.5, facecolor="#0f172a", edgecolor=SPINE_COL,
                 labelcolor="#e2e8f0")

# ── Panel 5: Mean trend per minute ───────────────────────────────────────────
ax5 = fig.add_subplot(gs[2, 1])
ax5.set_facecolor("#1e293b")
style_ax(ax5, "Panel E — Per-minute Mean  ±1 StdDev")
ax5.set_ylim(0.3, 0.7)

if minute_rows:
    mins  = [int(r["minute"])  for r in minute_rows]
    means = [float(r["mean"])  for r in minute_rows]
    sds   = [float(r["stddev"]) for r in minute_rows]
    means = np.array(means); sds = np.array(sds)

    ax5.plot(mins, means, color=C_OVERALL, linewidth=1.5, label="Mean/min")
    ax5.fill_between(mins, means-sds, means+sds, alpha=0.2, color=C_OVERALL,
                     label="±1 σ band")
    ax5.axhline(GAUSS_MU, color="#fbbf24", linestyle="--", linewidth=1.0,
                label="μ = 0.5")
    ax5.set_xlabel("Minute", color=TICK_COL, fontsize=8)
    ax5.set_ylabel("Mean value", **LABEL_KW)
    leg = ax5.legend(fontsize=7.5, facecolor="#0f172a", edgecolor=SPINE_COL,
                     labelcolor="#e2e8f0")

patches = [
    mpatches.Patch(color=C_OVERALL, label="60-min (Overall)"),
    mpatches.Patch(color=C_10MIN,   label="10-min blocks"),
    mpatches.Patch(color=C_1MIN,    label="1-min intervals"),
    mpatches.Patch(color="#fbbf24", label="Median (yellow line)"),
    mpatches.Patch(color="#f472b6", label="Mean (diamond marker)"),
    mpatches.Patch(color="#f97316", label="Tukey fences"),
]
fig.legend(handles=patches, loc="upper center", ncol=6,
           facecolor="#1e293b", edgecolor=SPINE_COL, labelcolor="#e2e8f0",
           fontsize=8, bbox_to_anchor=(0.5, 0.975))

fig.suptitle(
    "Gaussian Stream Statistics  —  μ=0.50, σ=0.15, 100 000 values/sec, 1 hour\n"
    "Central Tendency Analysis: Overall · 10-min · 1-min · PDF · Mean Trend",
    color="#f8fafc", fontsize=12, fontweight="bold", y=0.995)

plt.savefig("gaussian_boxplots.png", dpi=160, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print("Saved → gaussian_boxplots.png")