"""Main figures (PDF + 300-dpi PNG) in results/figures.

Fig 1  concept (A) and claim flow (B)
Fig 2  theta: OR of diagnosed prostate cancer per unit ln PSA by instrument set and cohort (A); analytic theta against
       screening intensity with the empirical cis-KLK estimate (B)
Fig 3  replicated claims: detection fraction with 95% interval (A); re-analysed and detection-corrected OR (B)
Fig 4  original grade vs new grade
Supplementary: S1 cohort-specific estimates for replicated claims; S2 observed vs PSA-predicted log OR (all claims)
Colours: orange = detection-related (filled = largely, open = partly), blue = not detection-related, grey = not replicated.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
ORANGE, BLUE, GREY = "#eb6834", "#2a78d6", "#8f8e8a"
GCOL = {"largely detection-related": ORANGE, "partly detection-related": ORANGE, "not detection-related": BLUE,
        "not replicated": GREY}
OUT_LAB = {"PCa_Wang2023": "Wang 2023 meta-analysis", "PCa_PRACTICAL": "PRACTICAL", "PCa_UKB": "UK Biobank",
           "PCa_FinnGen": "FinnGen R5"}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    if not name.startswith("FigS"):  # journal files: TIFF, 1200 dpi (line art), RGB, LZW-compressed
        import io
        from PIL import Image
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=1200, bbox_inches="tight")
        buf.seek(0)
        with Image.open(buf) as im:
            im.convert("RGB").save(os.path.join(FIG, f"{name}.tif"), compression="tiff_lzw", dpi=(1200, 1200))
    plt.close(fig)


def marker(ax, x, y, grade, size=5, **kw):
    c = GCOL[grade]
    face = "white" if grade == "partly detection-related" else c
    ax.plot(x, y, "o", ms=size, mfc=face, mec=c, mew=1.2, zorder=3, **kw)


# ---------------------------------------------------------------- Fig 1
def fig1(n):
    fig = plt.figure(figsize=(7.2, 3.0))
    ax = fig.add_axes([0, 0, 0.62, 1]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    box = dict(boxstyle="round,pad=0.35,rounding_size=0.15", fc="white", ec=INK, lw=0.8)

    def node(x, y, t, **kw):
        ax.text(x, y, t, ha="center", va="center", bbox={**box, **kw}, fontsize=8)

    node(0.8, 3.0, "Genetic\ninstruments")
    node(3.0, 3.0, "Exposure X")
    node(5.6, 4.9, "Latent prostate\ncancer")
    node(5.6, 1.1, "Serum PSA")
    node(8.2, 3.0, "Diagnosed prostate\ncancer (GWAS case)", ec=INK)
    arr = dict(arrowstyle="-|>", lw=0.9, color=INK, shrinkA=14, shrinkB=14, mutation_scale=9)

    def edge(a, b, t, tx, ty, color=INK, ls="-"):
        ax.annotate("", xy=b, xytext=a, arrowprops={**arr, "color": color, "ls": ls})
        ax.text(tx, ty, t, ha="center", va="center", fontsize=7.5, color=color)

    ax.annotate("", xy=(2.15, 3.0), xytext=(1.75, 3.0), arrowprops={**arr, "shrinkA": 0, "shrinkB": 0})
    edge((3.0, 3.0), (5.6, 4.9), "causal effect\n(target)", 3.75, 4.45, color=MUTED, ls="--")
    edge((3.0, 3.0), (5.6, 1.1), r"$\beta_{X\rightarrow PSA}$", 3.95, 1.75, color=ORANGE)
    edge((5.6, 4.9), (8.2, 3.0), "", 0, 0)
    edge((5.6, 1.1), (8.2, 3.0), r"$\theta$: PSA test, biopsy" + "\n" + "(detection)", 7.55, 1.45, color=ORANGE)
    ax.annotate("", xy=(5.6, 1.1), xytext=(5.6, 4.9), arrowprops={**arr, "color": MUTED, "lw": 0.6})
    ax.text(5.85, 3.0, "cancer\nraises PSA", fontsize=6.5, color=MUTED, va="center")
    ax.text(5.0, 0.05, r"Detection-driven log OR:  $E = \beta_{X\rightarrow PSA}\times\theta$;   "
            r"detection fraction $DF = E / \beta_{X\rightarrow PCa}$", ha="center", fontsize=7.5, color=INK)
    ax.text(0.05, 5.85, "A", fontsize=11, fontweight="bold", va="top")

    bx = fig.add_axes([0.64, 0.02, 0.36, 0.96]); bx.axis("off"); bx.set_xlim(0, 10); bx.set_ylim(0, 10)
    bx.text(0.0, 9.8, "B", fontsize=11, fontweight="bold", va="top")
    steps = [(f"{n['all']} MR associations with prostate cancer\n(Cui et al., PLoS Med 2024; search to Sept 2023)", 8.6),
             (f"{n['sig']} reported as significant", 6.4),
             (f"{n['re']} reanalyzed in European-ancestry GWAS", 4.2),
             (f"{n['rep']} replicated in Wang 2023\n(P < 0.05, same direction)", 2.0)]
    for t, y in steps:
        bx.text(4.6, y, t, ha="center", va="center", fontsize=7.5, bbox=box)
    for (_, y1), (_, y2) in zip(steps[:-1], steps[1:]):
        bx.annotate("", xy=(4.6, y2 + 0.55), xytext=(4.6, y1 - 0.55), arrowprops={**arr, "shrinkA": 0, "shrinkB": 0})
    bx.text(4.9, 5.3, f"{n['sig'] - n['re']} not reanalyzable:\nno summary statistics (4),\nEast Asian outcome (1)",
            ha="left", va="center", fontsize=6.5, color=MUTED)
    save(fig, "Fig1_concept_flow")


# ---------------------------------------------------------------- Fig 2
def fig2():
    th = pd.read_csv(os.path.join(TAB, "theta.tsv"), sep="\t")
    sim = pd.read_csv(os.path.join(TAB, "sim_theta.tsv"), sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.32))
    rows = []
    for ins, lab, meth in (("PSA_KLK", "cis-KLK", "IVW (LD-correlated)"), ("PSA_gw", "genome-wide", "IVW")):
        for o in OUT_LAB:
            r = th[(th.instruments == ins) & (th.outcome == o) & (th.method == meth)]
            if len(r):
                rows.append((f"{OUT_LAB[o]}", lab, r.iloc[0]))
    y = np.arange(len(rows))[::-1]
    for yi, (o, lab, r) in zip(y, rows):
        c = ORANGE if lab == "cis-KLK" else MUTED
        a.plot([r.lo, r.hi], [yi, yi], color=c, lw=1.2)
        a.plot(r.or_per_lnpsa, yi, "s" if lab == "cis-KLK" else "o", color=c, ms=5)
        a.text(0.99, yi, f"{r.or_per_lnpsa:.2f} ({r.lo:.2f}-{r.hi:.2f})", va="center", ha="right", fontsize=6.5,
               color=INK, transform=a.get_yaxis_transform())
    a.set_yticks(y, [f"{o}" for o, lab, r in rows], fontsize=7)
    a.set_xscale("log"); a.set_xlim(0.7, 18); a.axvline(1, color=GRID, lw=0.8, zorder=0)
    a.set_xticks([1, 2, 4], ["1", "2", "4"])
    a.set_xlabel("OR of diagnosed prostate cancer per unit ln PSA")
    a.text(-0.55, 1.02, "A", transform=a.transAxes, fontsize=11, fontweight="bold")
    nk = sum(1 for _, lab, _ in rows if lab == "cis-KLK")
    if len(rows) > nk:
        a.axhline(y[nk] + 0.5, color=GRID, lw=0.6)
    a.plot([], [], "s", color=ORANGE, ms=5, label="cis-KLK instruments")
    a.plot([], [], "o", color=MUTED, ms=5, label="genome-wide PSA instruments")
    a.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, handletextpad=0.3)
    a.set_ylim(-0.6, len(rows) - 0.4)
    base = sim[sim.scenario == "base"]
    sens = sim[sim.scenario != "base"]
    for T, ls in ((4.0, "-"), (3.0, "--")):
        s = base[base["T"] == T]
        b.plot(s.s, s.theta, color=INK, lw=1.4, ls=ls, label=f"biopsy threshold {T:g} ng/ml")
    for s0 in (0.3, 0.7):
        v = sens[(sens.s == s0) & (sens["T"] == 4.0)].theta
        b.plot([s0, s0], [v.min(), v.max()], color=MUTED, lw=3, alpha=0.35, solid_capstyle="butt")
    k = th[(th.instruments == "PSA_KLK") & (th.outcome == "PCa_Wang2023") & (th.method == "IVW (LD-correlated)")].iloc[0]
    b.axhspan(k.theta - 1.96 * k.se, k.theta + 1.96 * k.se, color=ORANGE, alpha=0.18, lw=0)
    b.axhline(k.theta, color=ORANGE, lw=1.2)
    b.text(0.98, k.theta - 2.2 * k.se, "cis-KLK θ, Wang 2023 (95% CI)", color=ORANGE, fontsize=7, ha="right",
           va="top")
    b.set_xlabel("Screening intensity (probability of a PSA test)")
    b.set_ylabel(r"$\theta$ (log OR per unit ln PSA)")
    b.set_xlim(0, 1); b.set_ylim(0, 1.45)
    b.legend(frameon=False, fontsize=7, loc="lower right")
    b.text(-0.2, 1.02, "B", transform=b.transAxes, fontsize=11, fontweight="bold")
    save(fig, "Fig2_theta")


# ---------------------------------------------------------------- Fig 3
def fig3(g):
    r = g[g.replicated].copy().sort_values("DF")
    n = len(r)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 0.28 * n + 1.1), sharey=True,
                               gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.08))
    y = np.arange(n)
    for yi, (_, x) in zip(y, r.iterrows()):
        lo, hi = max(x.DF_lo, -0.6), min(x.DF_hi, 1.6)
        a.plot([lo, hi], [yi, yi], color=GCOL[x.new_grade], lw=1.1)
        marker(a, x.DF, yi, x.new_grade)
        o1, l1, h1 = np.exp([x.b_pca, x.b_pca - 1.96 * x.se_pca, x.b_pca + 1.96 * x.se_pca])
        o2, l2, h2 = np.exp([x.b_corr, x.b_corr - 1.96 * x.se_corr, x.b_corr + 1.96 * x.se_corr])
        b.plot([l1, h1], [yi + 0.17, yi + 0.17], color=INK, lw=1.0)
        b.plot(o1, yi + 0.17, "o", color=INK, ms=4)
        b.plot([l2, h2], [yi - 0.17, yi - 0.17], color=GCOL[x.new_grade], lw=1.0)
        b.plot(o2, yi - 0.17, "D", mfc="white", mec=GCOL[x.new_grade], ms=4, mew=1.1)
    a.set_yticks(y, [f"{nice(x.exposure)} ({x.grade})" for _, x in r.iterrows()], fontsize=7)
    a.axvline(0, color=GRID, lw=0.8, zorder=0); a.axvline(0.5, color=GRID, lw=0.8, ls=":", zorder=0)
    a.set_xlim(-0.6, 1.6); a.set_xlabel("Detection fraction (95% interval, truncated at −0.6 and 1.6)")
    a.text(-0.02, 1.0 + 0.4 / n, "A", transform=a.transAxes, fontsize=11, fontweight="bold", ha="right")
    b.set_xscale("log"); b.axvline(1, color=GRID, lw=0.8, zorder=0)
    b.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.8, 0.9, 1.0, 1.1, 1.25, 1.5]))
    b.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    b.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_xlabel("OR per exposure unit (Wang 2023)")
    b.plot([], [], "o", color=INK, ms=4, label="reanalyzed")
    b.plot([], [], "D", mfc="white", mec=MUTED, ms=4, label="detection-corrected")
    b.legend(frameon=False, fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
    b.text(0.0, 1.0 + 0.4 / n, "B", transform=b.transAxes, fontsize=11, fontweight="bold")
    h = [plt.Line2D([], [], marker="o", ls="", mfc=ORANGE, mec=ORANGE, label="largely detection-related (DF ≥ 0.5)"),
         plt.Line2D([], [], marker="o", ls="", mfc="white", mec=ORANGE, mew=1.2, label="partly detection-related"),
         plt.Line2D([], [], marker="o", ls="", mfc=BLUE, mec=BLUE, label="not detection-related")]
    a.legend(handles=h, frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=1)
    a.set_ylim(-0.6, n - 0.4)
    save(fig, "Fig3_detection_fraction")


# ---------------------------------------------------------------- Fig 4
def fig4(g):
    order_old = ["Robust", "Probable", "Suggestive", "Insufficient"]
    order_new = ["largely detection-related", "partly detection-related", "not detection-related", "not replicated"]
    ct = pd.crosstab(g.grade, g.new_grade).reindex(index=order_old, columns=order_new, fill_value=0)
    fig, ax = plt.subplots(figsize=(5.2, 2.2))
    left = np.zeros(len(order_old))
    y = np.arange(len(order_old))[::-1]
    for k in order_new:
        v = ct[k].values
        face = "white" if k == "partly detection-related" else GCOL[k]
        hatch = "////" if k == "partly detection-related" else None
        ax.barh(y, v, left=left, height=0.62, color=face, edgecolor=GCOL[k] if face == "white" else "white",
                lw=1.0 if face == "white" else 1.5, hatch=hatch, label=k)
        for yi, l, vv in zip(y, left, v):
            if vv > 0:
                ax.text(l + vv / 2, yi, str(vv), ha="center", va="center", fontsize=7,
                        color="white" if face not in ("white",) and k != "not replicated" else INK)
        left += v
    ax.set_yticks(y, [f"{o} (n={int(ct.loc[o].sum())})" for o in order_old])
    ax.set_xlabel("Number of reanalyzed claims")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Original grade (Cui 2024) and reanalysis result", fontsize=8, loc="left")
    save(fig, "Fig4_regrade")


# ---------------------------------------------------------------- Supplementary
def figS1(d, g):
    rep = g[g.replicated].sort_values("DF")
    keys = list(rep.claim)
    n = len(keys)
    fig, ax = plt.subplots(figsize=(4.8, 0.36 * n + 0.9))
    offs = {"PCa_Wang2023": 0.27, "PCa_PRACTICAL": 0.09, "PCa_UKB": -0.09, "PCa_FinnGen": -0.27}
    cols = {"PCa_Wang2023": INK, "PCa_PRACTICAL": "#2a78d6", "PCa_UKB": "#eb6834", "PCa_FinnGen": "#1baf7a"}
    for yi, c in enumerate(keys):
        for o, off in offs.items():
            r = d[(d.claim == c) & (d.outcome == o)]
            if r.empty:
                continue
            r = r.iloc[0]
            z = r.b_pca / abs(rep[rep.claim == c].b_pca.iloc[0])
            s = r.se_pca / abs(rep[rep.claim == c].b_pca.iloc[0])
            ax.plot([z - 1.96 * s, z + 1.96 * s], [yi + off] * 2, color=cols[o], lw=0.9)
            ax.plot(z, yi + off, "o", color=cols[o], ms=3.2)
    ax.set_yticks(range(n), [nice(rep[rep.claim == c].exposure.iloc[0]) for c in keys], fontsize=7)
    ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    ax.set_xlim(-4, 4)
    ax.set_xlabel("Estimate relative to |Wang 2023 estimate| (95% CI, clipped)")
    for o in offs:
        ax.plot([], [], "o", color=cols[o], ms=4, label=OUT_LAB[o])
    ax.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    save(fig, "FigS1_cohorts")


def figS2(g):
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    z_obs = g.b_pca / g.se_pca
    for _, x in g.iterrows():
        sc = abs(x.b_pca) if x.replicated else np.nan
    # standardise each claim by the SE of its observed estimate so that units cancel
    xs, ys = g.E / g.se_pca, g.b_pca / g.se_pca
    for gr in ["not replicated", "not detection-related", "partly detection-related", "largely detection-related"]:
        m = g.new_grade == gr
        face = "white" if gr == "partly detection-related" else GCOL[gr]
        ax.scatter(xs[m], ys[m], s=22, facecolor=face, edgecolor=GCOL[gr], lw=1.0, label=gr, zorder=3)
    lim = max(abs(xs).max(), abs(ys).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=0.7, ls="--")
    ax.axhline(0, color=GRID, lw=0.6); ax.axvline(0, color=GRID, lw=0.6)
    ax.set_xlabel("PSA-predicted log OR / SE(observed)")
    ax.set_ylabel("Observed log OR / SE(observed)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    save(fig, "FigS2_observed_vs_predicted")


if __name__ == "__main__":
    g = pd.read_csv(os.path.join(TAB, "regrade_PSA_KLK.tsv"), sep="\t")
    d = pd.read_csv(os.path.join(TAB, "detection_PSA_KLK.tsv"), sep="\t")
    cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")
    n = dict(all=145, sig=len(cm), re=int((cm.inst_mode != "exclude").sum()), rep=int(g.replicated.sum()))
    fig1(n); fig2(); fig3(g); fig4(g); figS1(d, g); figS2(g)
    print("figures written to", FIG)
