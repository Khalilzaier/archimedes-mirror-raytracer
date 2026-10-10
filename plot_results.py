"""Trace les figures de l'article à partir des CSV produits par ./archimede.

Usage (depuis la racine du dépôt) :  python3 scripts/plot_results.py
Les figures sont écrites dans figures/.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RES = Path("results")
FIG = Path("figures")
FIG.mkdir(exist_ok=True)

# Erreur de visée = grandeur ordonnée : une seule teinte, du clair au foncé.
BLEUS = ["#8fb8ea", "#4f8fde", "#2a62b0", "#0f3466"]
ENCRE = "#1f1f1d"
ENCRE_2 = "#5a5955"
GRILLE = "#e4e3df"
ROUGE_SEUIL = "#b3261e"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": ENCRE_2,
    "axes.labelcolor": ENCRE,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRILLE,
    "grid.linewidth": 0.8,
    "xtick.color": ENCRE_2,
    "ytick.color": ENCRE_2,
    "lines.linewidth": 2.0,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


def lire(nom):
    return np.genfromtxt(RES / nom, delimiter=",", names=True)


# ---------------------------------------------------------------------------
# 1. Validation : un miroir parfait contre la théorie
# ---------------------------------------------------------------------------
def fraction_disque(hx, hy, alpha, n=1200):
    """Fraction du disque solaire (rayon angulaire alpha) vue à travers un
    rectangle de demi-côtés angulaires hx, hy centré sur le soleil réfléchi."""
    t = np.linspace(-alpha, alpha, n)
    X, Y = np.meshgrid(t, t)
    dans_disque = X**2 + Y**2 <= alpha**2
    dans_rect = (np.abs(X) <= hx) & (np.abs(Y) <= hy)
    return (dans_disque & dans_rect).sum() / dans_disque.sum()


def figure_validation():
    v = lire("validation.csv")
    D = v["D"]
    a, alpha, rho, dni = v["cote"][0], v["alpha"][0], v["rho"][0], v["dni"][0]
    # l'angle d'incidence dépend (un peu) de D : la cible est 1 m plus bas que le miroir
    Dth = np.linspace(D.min(), D.max(), 300)
    theta_th = np.interp(Dth, D, v["theta_inc_rad"])
    L = np.sqrt(Dth**2 + 1.0)                       # miroir 1 m plus haut que la cible
    cosg = Dth / L
    F = np.array([fraction_disque(a / (2 * l), a * np.cos(t) / (2 * l), alpha)
                  for l, t in zip(L, theta_th)])
    E_th = rho * dni * cosg * F

    # comparaison point par point
    Lp = np.sqrt(D**2 + 1.0)
    Fp = np.array([fraction_disque(a / (2 * l), a * np.cos(t) / (2 * l), alpha)
                   for l, t in zip(Lp, v["theta_inc_rad"])])
    Ep = rho * dni * (D / Lp) * Fp
    ecart = np.abs(v["flux_centre"] - Ep) / Ep * 100

    theta_loin = v["theta_inc_rad"][-1]
    Dstar = a * np.cos(theta_loin) / (2 * alpha)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(Dth, E_th, color=BLEUS[3], label="Théorie (recouvrement disque solaire / miroir)")
    ax.plot(D, v["flux_centre"], "o", ms=6, mfc="white", mec=BLEUS[1], mew=1.8,
            label="Lancer de rayons (10 millions de rayons)")
    ax.axvline(Dstar, color=ENCRE_2, lw=1, ls="--")
    ax.text(Dstar + 4, 120, f"$D^* = a\\,\\cos\\theta\\,/\\,2\\alpha \\approx {Dstar:.0f}$ m",
            color=ENCRE_2, fontsize=9)
    ax.set_xlabel("Distance miroir - coque (m)")
    ax.set_ylabel("Flux au centre de la tache (W/m²)")
    ax.set_title("Un miroir plan seul : validation du code", loc="left", color=ENCRE)
    ax.set_ylim(0, 620)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    fig.savefig(FIG / "fig1_validation.png")
    plt.close(fig)
    for Di, e in zip(D, ecart):
        print(f"  D = {Di:5.0f} m  ecart = {e:.2f} %")
    print(f"Validation : ecart max {ecart.max():.2f} %, moyen {ecart.mean():.2f} %, D* = {Dstar:.1f} m")


# ---------------------------------------------------------------------------
# 2. Carte de flux sur la coque
# ---------------------------------------------------------------------------
def figure_carte():
    chemin = RES / "flux_map.csv"
    entete = chemin.open().readline().lstrip("# ").split()
    p = dict(kv.split("=") for kv in entete)
    demi, pas = float(p["demi_largeur"]), float(p["pas"])
    M = np.loadtxt(chemin, delimiter=",", comments="#") / 1000.0  # kW/m2
    # lissage 3x3 (comme dans le code C++) pour les isolignes
    K = np.pad(M, 1)
    L = sum(K[1 + dj:K.shape[0] - 1 + dj, 1 + di:K.shape[1] - 1 + di]
            for dj in (-1, 0, 1) for di in (-1, 0, 1)) / 9.0

    zoom = 2.0
    ext = [-demi, demi, -demi, demi]
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    im = ax.imshow(M, origin="lower", extent=ext, cmap="inferno", vmin=0,
                   vmax=np.percentile(M, 99.9), interpolation="nearest")
    xs = np.linspace(-demi + pas / 2, demi - pas / 2, M.shape[0])
    cs = ax.contour(xs, xs, L, levels=[10, 20], colors=["#cfe3ff", "white"],
                    linewidths=[1.2, 1.8], linestyles=["--", "-"])
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color="#cfe3ff", ls="--", lw=1.2, label="10 kW/m² (incandescence)"),
                       Line2D([], [], color="white", lw=1.8, label="20 kW/m² (inflammation)")],
              loc="lower left", fontsize=8, frameon=True, facecolor="#1f1f1d", edgecolor="none",
              labelcolor="white")
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(-zoom, zoom)
    ax.set_xlabel("Position le long de la coque (m)")
    ax.set_ylabel("Hauteur relative au point visé (m)")
    ax.grid(False)
    ax.set_title(f"{p['N']} miroirs à {float(p['D']):.0f} m, visée ±{float(p['aim_deg'])}°",
                 loc="left", color=ENCRE)
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Flux incident (kW/m²)")
    fig.savefig(FIG / "fig2_carte_flux.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3. Nombre de miroirs requis selon la distance et l'erreur de visée
# ---------------------------------------------------------------------------
def figure_requis():
    r = lire("sweep_required.csv")
    visees = np.unique(r["aim_deg"])
    nmax = int(lire_nmax())
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for k, a in enumerate(visees):
        sel = r["aim_deg"] == a
        Ds = np.unique(r["D"][sel])
        med = []
        for D in Ds:
            vals = r["N_req_flamme"][sel & (r["D"] == D)]
            vals = np.where(vals < 0, np.nan, vals)
            med.append(np.nan if np.isnan(vals).any() else np.median(vals))
        med = np.array(med)
        ax.plot(Ds, med, "-o", ms=4, color=BLEUS[k])
        ok = ~np.isnan(med)
        if ok.any():
            i = np.where(ok)[0][-1]
            ax.annotate(f"±{a:g}°", (Ds[i], med[i]), xytext=(6, 0), textcoords="offset points",
                        color=ENCRE, fontsize=9, va="center")
    ax.axhline(nmax, color=ENCRE_2, lw=1, ls=":")
    ax.text(152, nmax * 1.05, f"limite simulée : {nmax} miroirs", color=ENCRE_2, fontsize=8, ha="right")
    ax.set_yscale("log")
    ticks = [40, 60, 100, 200, 300, 500]
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(33, nmax * 1.25)
    ax.set_xlabel("Distance du navire (m)")
    ax.set_ylabel("Miroirs nécessaires (médiane)")
    ax.set_title("Miroirs requis pour 20 kW/m² selon l'erreur de visée", loc="left", color=ENCRE)
    fig.savefig(FIG / "fig3_miroirs_requis.png")
    plt.close(fig)


def lire_nmax():
    p = lire("sweep_peak.csv")
    return p["N"].max()


# ---------------------------------------------------------------------------
# 4. Pic de flux selon le nombre de miroirs (D fixé)
# ---------------------------------------------------------------------------
def figure_pic(D_choisi=50):
    p = lire("sweep_peak.csv")
    visees = np.unique(p["aim_deg"])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for k, a in enumerate(visees):
        sel = (p["aim_deg"] == a) & (p["D"] == D_choisi)
        Ns = np.unique(p["N"][sel])
        moy = np.array([p["peak"][sel & (p["N"] == n)].mean() for n in Ns]) / 1000
        ax.plot(Ns, moy, color=BLEUS[k])
        ax.annotate(f"±{a:g}°", (Ns[-1], moy[-1]), xytext=(4, 0), textcoords="offset points",
                    color=ENCRE, fontsize=9, va="center")
    ax.axhline(20, color=ROUGE_SEUIL, lw=1.2)
    ax.text(5, 21.5, "auto-inflammation du bois (20 kW/m²)", color=ROUGE_SEUIL, fontsize=8)
    ax.axhline(10, color=ROUGE_SEUIL, lw=1, ls="--")
    ax.text(5, 11.5, "incandescence (10 kW/m²)", color=ROUGE_SEUIL, fontsize=8)
    ax.set_xlabel("Nombre de miroirs")
    ax.set_ylabel("Pic de flux sur la coque (kW/m²)")
    ax.set_title(f"Navire à {D_choisi} m : le flux plafonne quand la visée est imprécise",
                 loc="left", color=ENCRE)
    fig.savefig(FIG / "fig4_pic_vs_N.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5. Effet de l'heure de la journée
# ---------------------------------------------------------------------------
def figure_journee():
    d = lire("day.csv")
    hs = np.unique(d["heure"])
    med, theta = [], []
    for h in hs:
        sel = d["heure"] == h
        vals = d["N_req_flamme"][sel]
        med.append(np.nan if (vals < 0).any() else np.median(vals))
        theta.append(d["theta_moyen_deg"][sel].mean())
    med = np.array(med)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(hs, med, "-o", ms=4, color=BLEUS[2])
    imp = np.isnan(med)
    if imp.any():
        for h in hs[imp]:
            ax.axvspan(h - 0.25, h + 0.25, color=GRILLE, lw=0)
        xm = hs[imp].mean()
        ax.text(xm, np.nanmin(med) + 0.5 * (np.nanmax(med) - np.nanmin(med)),
                "seuil jamais atteint\n(moins de 500 miroirs)", ha="center", va="center",
                rotation=90, color=ENCRE_2, fontsize=8)
    ax.set_xlabel("Heure solaire (Syracuse, début juillet)")
    ax.set_ylabel("Miroirs nécessaires (médiane)")
    ax.text(11, np.nanmin(med) + 0.6 * (np.nanmax(med) - np.nanmin(med)), "navire à 50 m, visée ±0.25°", color=ENCRE_2, fontsize=9, ha="center")
    ax.set_title("Fenêtre favorable en fin de matinée, impossible le soir", loc="left", color=ENCRE)
    fig.savefig(FIG / "fig5_heure.png")
    plt.close(fig)


if __name__ == "__main__":
    for nom, f in [("validation.csv", figure_validation), ("flux_map.csv", figure_carte),
                   ("sweep_required.csv", figure_requis), ("sweep_peak.csv", figure_pic),
                   ("day.csv", figure_journee)]:
        if (RES / nom).exists():
            f()
    print("Figures -> figures/")
