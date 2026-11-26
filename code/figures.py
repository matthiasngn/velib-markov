"""Figures de l'article. Usage : python3 -I figures.py <dossier_sorties> <dossier_figures>"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SORTIES, FIG = sys.argv[1], sys.argv[2]
S = np.load(f"{SORTIES}/series.npz")
R = json.load(open(f"{SORTIES}/resultats.json"))

ENCRE, GRIS, CLAIR = "#1f2328", "#6b7280", "#c9ced6"
BLEU, ORANGE = "#2a78d6", "#eb6834"
COL, LARGE = 3.31, 6.85                 # largeurs en pouces (colonne, page)

plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.5, "axes.edgecolor": GRIS, "axes.labelcolor": ENCRE,
    "xtick.color": GRIS, "ytick.color": GRIS, "xtick.labelcolor": ENCRE, "ytick.labelcolor": ENCRE,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5, "xtick.minor.width": 0.4,
    "axes.grid": True, "grid.color": "#e6e8eb", "grid.linewidth": 0.4, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.handlelength": 1.6, "legend.borderaxespad": 0.2, "legend.labelspacing": 0.25,
    "lines.linewidth": 0.9, "savefig.bbox": "tight", "savefig.pad_inches": 0.015, "pdf.fonttype": 42,
    "axes.titlesize": 8, "axes.titlelocation": "left", "axes.titlepad": 3,
})


def virgule(x, nd=1):
    return f"{x:.{nd}f}".replace(".", "{,}")


def bande(ax, t, q05, q95, moy, label_bande="bande à 90 %", label_moy="espérance"):
    ax.fill_between(t, q05, q95, step="post", color=BLEU, alpha=0.16, lw=0, label=label_bande)
    ax.plot(t, moy, color=BLEU, lw=0.9, label=label_moy)


def heures_jour(ax, debut=0, fin=24):
    ax.set_xlim(debut, fin)
    ax.set_xticks(range(debut, fin + 1, 6))
    lab = ax.set_xticklabels([f"{h % 24} h" for h in range(debut, fin + 1, 6)])
    lab[0].set_ha("left"); lab[-1].set_ha("right")


# ---------------------------------------------------------------- 1. données
fig, axs = plt.subplots(1, 2, figsize=(LARGE, 1.75), gridspec_kw=dict(width_ratios=[2, 1], wspace=0.12))
jours = ["mar. 2", "mer. 3", "jeu. 4", "ven. 5", "sam. 6"]
for ax, t, x, cap, nbj, titre in [
        (axs[0], S["gdl_t"], S["gdl_x"], 45, 4, "(a) Gare de Lyon, $C=45$"),
        (axs[1], S["nat_t"], S["nat_x"], 36, 2, "(b) Nation, $C=36$")]:
    ok = (t >= 0) & (t <= 24 * nbj)
    ax.step(t[ok], x[ok], where="post", color=ENCRE, lw=0.6)
    ax.axhline(cap, color=GRIS, lw=0.6, ls=(0, (3, 2)))
    ax.set_xlim(0, 24 * nbj)
    ax.set_xticks(range(0, 24 * nbj + 1, 24))
    ax.set_xticklabels(jours[:nbj + 1])
    ax.set_xticks(range(0, 24 * nbj + 1, 6), minor=True)
    ax.set_ylim(0, 50)
    ax.set_title(titre)
axs[0].set_ylabel("vélos disponibles")
axs[0].text(1, 45.9, "capacité", color=GRIS, fontsize=6.5, va="bottom")
fig.savefig(f"{FIG}/fig_donnees.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 2. stationnaire
st = R["stationnaire"]["gdl"]
fig, axs = plt.subplots(2, 1, figsize=(COL, 3.05), gridspec_kw=dict(hspace=0.55))
ax = axs[0]
n = np.arange(46)
ax.bar(n - 0.21, S["st_pi"], width=0.42, color=ORANGE, lw=0, label=r"loi stationnaire $\pi$")
ax.bar(n + 0.21, S["st_emp"], width=0.42, color=ENCRE, lw=0, label="occupation observée")
ax.set_xlim(-1, 46)
ax.set_xlabel("nombre de vélos $n$")
ax.set_ylabel("proportion")
ax.set_ylim(0, 0.36)
ax.legend(loc="upper center", ncol=2, columnspacing=1.0)
ax.set_title("(a) loi stationnaire et occupation observée")
ax = axs[1]
tt = 9 + np.arange(len(S["st_x"])) / 60
bande(ax, tt, S["st_q05"], S["st_q95"], S["st_moy"], label_moy="espérance (régime transitoire)")
ax.step(tt, S["st_x"], where="post", color=ENCRE, lw=0.6, label="observé")
ax.axhline(st["moy_stat"], color=ORANGE, lw=0.9, ls=(0, (3, 2)), label=r"$\mathbb{E}_\pi[X]$")
ax.set_xlim(9, 19)
ax.set_xticks(range(9, 20, 2))
ax.set_xticklabels([f"{h} h" for h in range(9, 20, 2)])
ax.set_ylim(0, 46)
ax.set_ylabel("vélos")
ax.text(9.15, 38.5, "observé", fontsize=7, color=ENCRE)
ax.text(18.9, st["moy_stat"] + 1.2, r"$\mathbb{E}_\pi[X]$", fontsize=7, color=ENCRE, ha="right")
ax.text(11.2, 6.8, "espérance et bande à 90 % de la loi de $X_t$", fontsize=7, color=ENCRE)
ax.set_title("(b) trajectoire du mardi, 9 h – 19 h")
fig.savefig(f"{FIG}/fig_stationnaire.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 3. forçage
fig, ax = plt.subplots(figsize=(COL, 1.85))
h = np.arange(1, 24, 2)
for taux, ic, couleur, dec, lab in [(S["lam_all"], S["lam_ic"], BLEU, -0.3, r"dépôts $\hat\lambda_k$"),
                                     (S["mu_all"], S["mu_ic"], ORANGE, 0.3, r"prises $\hat\mu_k$")]:
    bas = np.where(np.isnan(ic[:, 1]), 0, taux - ic[:, 0])
    haut = np.where(np.isnan(ic[:, 1]), 0, ic[:, 1] - taux)
    ax.errorbar(h + dec, taux, yerr=[bas, haut], fmt="o", ms=2.6, color=couleur, elinewidth=0.7, capsize=0, label=lab)
ax.set_xlim(0, 24)
ax.set_xticks(range(0, 25, 4))
ax.set_xticklabels([f"{x} h" for x in range(0, 25, 4)])
ax.set_ylabel(r"taux (h$^{-1}$)")
ax.set_ylim(0, 13.5)
ax.legend(loc="upper right")
fig.savefig(f"{FIG}/fig_forcage.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 4. vitesse de convergence
fig, ax = plt.subplots(figsize=(COL, 1.9))
t = S["dist_t"]
for cle, couleur, lab in [("equil", BLEU, "8–10 h"), ("pointe", ORANGE, "6–8 h")]:
    d = R["dist"][cle]
    ax.semilogy(t, np.maximum(S[f"dist_{cle}"], 1e-12), color=couleur, lw=0.9,
                label=f"{lab} ($\\hat\\rho={virgule(d['lam'] / d['mu'], 2 if cle == 'equil' else 1)}$)")
    ax.semilogy(t, S[f"borne_{cle}"], color=couleur, lw=0.6, ls=(0, (3, 2)))
ax.axhline(0.25, color=GRIS, lw=0.5, ls=":")
ax.text(47.5, 0.29, r"$\varepsilon=1/4$", ha="right", va="bottom", fontsize=7, color=GRIS)
ax.set_ylim(1e-4, 3)
ax.set_xlim(0, 48)
ax.set_xticks(range(0, 49, 12))
ax.set_xlabel("$t$ (h)")
ax.set_ylabel(r"$\|\delta_0 e^{tQ}-\pi\|_{\mathrm{VT}}$")
ax.plot([], [], color=GRIS, lw=0.6, ls=(0, (3, 2)), label="borne supérieure (th. 3.4)")
ax.legend(loc="lower right")
fig.savefig(f"{FIG}/fig_melange.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 5. prévision à 2 h
fig, axs = plt.subplots(1, 2, figsize=(COL, 1.75), sharey=True, gridspec_kw=dict(wspace=0.12))
for ax, cle, h0, titre in [(axs[0], "matin", 6, "(a) mardi 6–8 h"), (axs[1], "soir", 20, "(b) mardi 20–22 h")]:
    tt = h0 + np.arange(121) / 60
    for sim in S[f"c_{cle}_sims"]:
        ax.step(tt, sim, where="post", color=CLAIR, lw=0.4)
    bande(ax, tt, S[f"c_{cle}_q05"], S[f"c_{cle}_q95"], S[f"c_{cle}_moy"])
    ax.step(tt, S[f"c_{cle}_x"], where="post", color=ENCRE, lw=0.7, label="observé")
    ax.set_xlim(h0, h0 + 2)
    ax.set_xticks([h0, h0 + 1, h0 + 2])
    lab = ax.set_xticklabels([f"{h0} h", f"{h0 + 1} h", f"{h0 + 2} h"])
    lab[0].set_ha("left"); lab[-1].set_ha("right")
    ax.set_title(titre)
axs[0].set_ylim(0, 46)
axs[0].set_ylabel("vélos")
axs[0].plot([], [], color=CLAIR, lw=0.8, label="trajectoires simulées")
h_, l_ = axs[0].get_legend_handles_labels()
ordre = [3, 0, 1, 2]
fig.legend([h_[i] for i in ordre], [l_[i] for i in ordre], loc="lower center", ncol=2,
           bbox_to_anchor=(0.5, -0.2), columnspacing=1.0)
fig.savefig(f"{FIG}/fig_court.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 6. prévision à 24 h
fig, axs = plt.subplots(1, 4, figsize=(LARGE, 1.75), sharey=True, gridspec_kw=dict(wspace=0.08))
noms = {"02": "mardi 2", "03": "mercredi 3", "04": "jeudi 4", "05": "vendredi 5"}
for ax, j in zip(axs, ["02", "03", "04", "05"]):
    x = S[f"j_gdl_{j}_x"]
    tt = np.arange(len(x)) / 60
    bande(ax, tt, S[f"j_gdl_{j}_creneaux_q05"], S[f"j_gdl_{j}_creneaux_q95"], S[f"j_gdl_{j}_creneaux_moy"],
          label_moy="espérance (créneaux)")
    ax.plot(tt, S[f"j_gdl_{j}_constant_moy"], color=ORANGE, lw=0.8, ls=(0, (3, 2)), label="espérance (taux constants)")
    ax.step(tt, x, where="post", color=ENCRE, lw=0.6, label="observé")
    heures_jour(ax)
    ax.set_title(f"{noms[j]} juin")
axs[0].set_ylim(0, 46)
axs[0].set_ylabel("vélos")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.1))
fig.savefig(f"{FIG}/fig_jour.pdf")
plt.close(fig)

# ---------------------------------------------------------------- 7. régime périodique
fig, ax = plt.subplots(figsize=(LARGE, 1.7))
x = S["q4_x"]
tt = np.arange(len(x)) / 60
bande(ax, tt, S["q4_q05"], S["q4_q95"], S["q4_moy"], label_moy=r"espérance, $X_0=1$")
per = np.tile(S["per_moy"][:-1], 4)[:len(x)]
ax.plot(tt, per, color=ORANGE, lw=0.8, ls=(0, (3, 2)), label=r"régime périodique $\mathbb{E}_{\nu}[X_t]$")
ax.step(tt, x, where="post", color=ENCRE, lw=0.6, label="observé")
ax.set_xlim(0, 96)
ax.set_xticks(range(0, 97, 24))
ax.set_xticklabels(jours)
ax.set_xticks(range(0, 97, 6), minor=True)
ax.set_ylim(0, 46)
ax.set_ylabel("vélos")
ax.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.17))
fig.savefig(f"{FIG}/fig_periodique.pdf")
plt.close(fig)
print("figures OK")
