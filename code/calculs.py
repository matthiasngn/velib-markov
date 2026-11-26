"""Tous les calculs de l'article. Sorties : resultats.json (nombres) et series.npz (courbes).

Usage : python3 -I calculs.py <dossier_scripts> <dossier_donnees> <dossier_sortie>
"""
import csv
import json
import sys
from datetime import datetime, timedelta

import numpy as np
from scipy.linalg import expm
from scipy.stats import binomtest, chi2, poisson, skellam

sys.path.insert(0, sys.argv[1])
import donnees as D
import velib as V

DATA, OUT = sys.argv[2], sys.argv[3]
R = {}          # nombres
S = {}          # séries pour les figures
HC, NB = 2, 12
rng = np.random.default_rng(2026)


def vt(p, q):
    return 0.5 * float(np.abs(p - q).sum())


def bande(lois):
    etats = np.arange(lois.shape[1])
    moy = lois @ etats
    q05 = np.array([V.quantile(p, 0.05) for p in lois])
    q95 = np.array([V.quantile(p, 0.95) for p in lois])
    return moy, q05, q95


def comptes_creneaux(inst, velos, libres, jour):
    out = np.zeros((NB, 4))
    for k in range(NB):
        idx = D.indices(inst, jour + timedelta(hours=HC * k), jour + timedelta(hours=HC * (k + 1)))
        if len(idx) >= 2:
            out[k] = V.estimer_taux(velos[idx], libres[idx])
    return out


def taux_depuis(s):
    lam = np.divide(s[:, 0], s[:, 1], out=np.zeros(len(s)), where=s[:, 1] > 0)
    mu = np.divide(s[:, 2], s[:, 3], out=np.zeros(len(s)), where=s[:, 3] > 0)
    return lam, mu


def ic(nb, expo):
    if nb == 0:
        return 0.0, 0.0, float("nan")
    t, a, b = V.intervalle(nb, expo)
    return float(t), float(a), float(b)


def ic_rho(n_d, t_d, n_p, t_p):
    rho = (n_d / t_d) / (n_p / t_p)
    z = 1.96 * np.sqrt(1 / n_d + 1 / n_p)
    return float(rho), float(rho * np.exp(-z)), float(rho * np.exp(z))


def t_melange(C, lam, mu, eps=0.25, tmax=2000.0):
    Q = V.generateur(C, lam, mu)
    pi = V.loi_stationnaire(C, lam, mu)

    def d(t):
        return 0.5 * np.abs(expm(t * Q) - pi).sum(axis=1).max()

    a, b = 0.0, tmax
    assert d(b) <= eps
    for _ in range(50):
        m = (a + b) / 2
        if d(m) > eps:
            a = m
        else:
            b = m
    return float(b)


def rv(N, T):
    ok = T > 0
    N, T = N[ok], T[ok]
    tot = N.sum() / T.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        terme = np.where(N > 0, N * np.log((N / T) / tot), 0.0)
    stat = 2 * terme.sum()
    ddl = int(ok.sum() - 1)
    return float(stat), ddl, float(chi2.sf(stat, ddl))


def nb_lignes(nom):
    with open(f"{DATA}/{nom}") as f:
        return sum(1 for _ in csv.DictReader(f))


# =====================================================================
# 0. Vérifications des formules fermées
# =====================================================================
ecart = dict(spectre=0.0, vecteurs=0.0, atteinte=0.0, ruine=0.0)
for (C_, lam, mu) in [(45, 3.0, 3.3), (45, 10.05, 0.5), (36, 2.0, 2.0), (7, 1.0, 2.5), (12, 4.0, 4.0)]:
    Q = V.generateur(C_, lam, mu)
    rho = lam / mu
    n = np.arange(C_ + 1)
    theta = np.arange(1, C_ + 1) * np.pi / (C_ + 1)
    s_k = -(lam + mu) + 2 * np.sqrt(lam * mu) * np.cos(theta)
    d = rho ** (n / 2)
    Ssym = d[:, None] * Q / d[None, :]          # D^{1/2} Q D^{-1/2}, D = diag(pi)
    assert np.abs(Ssym - Ssym.T).max() < 1e-9 * np.abs(Q).max()
    num = np.linalg.eigvalsh((Ssym + Ssym.T) / 2)
    ecart["spectre"] = max(ecart["spectre"], np.abs(num - np.sort(np.concatenate([[0.0], s_k]))).max())
    for j in range(C_):
        g = np.sin((n + 1) * theta[j]) - rho ** (-0.5) * np.sin(n * theta[j])
        ecart["vecteurs"] = max(ecart["vecteurs"], np.abs(Ssym @ g - s_k[j] * g).max() / np.abs(g).max())
    assert abs(V.trou_spectral(C_, lam, mu) + s_k[0]) < 1e-12
    if max(rho, 1 / rho) ** C_ > 1e10:
        continue                                # systeme lineaire trop mal conditionne
    m_num = V.temps_moyen_atteinte(Q, [0])
    if abs(rho - 1) > 1e-12:
        m_th = (n - rho ** (C_ + 1 - n) * (1 - rho ** n) / (1 - rho)) / (mu * (1 - rho))
        h_th = (1 - rho ** (-n)) / (1 - rho ** (-C_))
    else:
        m_th = n * (2 * C_ - n + 1) / (2 * mu)
        h_th = n / C_
    ecart["atteinte"] = max(ecart["atteinte"], np.abs(m_num - m_th).max() / m_th.max())
    h_num = np.zeros(C_ + 1); h_num[C_] = 1.0
    h_num[1:C_] = np.linalg.solve(Q[1:C_, 1:C_], -Q[1:C_, C_])
    ecart["ruine"] = max(ecart["ruine"], np.abs(h_num - h_th).max())
R["verif"] = {k: float(v) for k, v in ecart.items()}
print("vérifications :", R["verif"])

# =====================================================================
# 1. Données
# =====================================================================
C = 45
inst, meca, elec, libres = D.semaine_gdl(DATA)
velos = meca + elec
JOURS = [datetime(2026, 6, j) for j in (2, 3, 4, 5)]
Cn = 36
inst_n, meca_n, elec_n, libres_n = D.semaine_nation(DATA)
velos_n = meca_n + elec_n
JN = [datetime(2026, 6, 2), datetime(2026, 6, 3)]
WE = {}
for st in ("gare_de_lyon", "nation", "parc_des_princes", "claude_terrasse"):
    i_w, m_w, e_w, l_w = D.weekend(DATA, st)
    WE[st] = (i_w, m_w + e_w, l_w)


def resume(nom, ii, vv, ll, cap, fichiers):
    d = np.diff(vv)
    R.setdefault("donnees", {})[nom] = dict(
        debut=ii[0].isoformat(), fin=ii[-1].isoformat(), heures=(len(ii) - 1) / 60,
        minutes=len(ii), capacite=cap,
        depots=int(np.maximum(d, 0).sum()), prises=int(np.maximum(-d, 0).sum()),
        mini=int(vv.min()), maxi=int(vv.max()),
        min_pleine=int((ll == 0).sum()), min_vide=int((vv == 0).sum()))


def nb_brut(fichiers):
    vus = set()
    for nom in fichiers:
        with open(f"{DATA}/{nom}") as f:
            for ligne in csv.DictReader(f):
                vus.add(ligne["instant"])
    return len(vus)


for nom, ii, vv, ll, cap, fichiers in [
    ("gdl_semaine", inst, velos, libres, 45, ["gare_de_lyon_20260601.csv", "gare_de_lyon_20260603.csv"]),
    ("nation_semaine", inst_n, velos_n, libres_n, 36, ["nation_20260601.csv"]),
    ("gdl_we", *WE["gare_de_lyon"], 45, ["gare_de_lyon_20260523.csv"]),
    ("nation_we", *WE["nation"], 36, ["nation_20260523.csv"]),
    ("princes_we", *WE["parc_des_princes"], 32, ["parc_des_princes_20260523.csv"]),
    ("terrasse_we", *WE["claude_terrasse"], 46, ["claude_terrasse_20260523.csv"]),
]:
    resume(nom, ii, vv, ll, cap, fichiers)
    R["donnees"][nom]["completees"] = len(ii) - nb_brut(fichiers)
    print(nom, R["donnees"][nom])

# comptages : toutes variations / unitaires seulement / par type de vélo
dv = np.diff(velos); dm = np.diff(meca); de = np.diff(elec)
R["comptage"] = dict(
    depots=int(np.maximum(dv, 0).sum()), prises=int(np.maximum(-dv, 0).sum()),
    depots_unit=int((dv == 1).sum()), prises_unit=int((dv == -1).sum()),
    depots_type=int((np.maximum(dm, 0) + np.maximum(de, 0)).sum()),
    prises_type=int((np.maximum(-dm, 0) + np.maximum(-de, 0)).sum()))
R["comptage"]["perte_unit_depots"] = 1 - R["comptage"]["depots_unit"] / R["comptage"]["depots"]
R["comptage"]["perte_unit_prises"] = 1 - R["comptage"]["prises_unit"] / R["comptage"]["prises"]
R["comptage"]["gain_type"] = (R["comptage"]["depots_type"] + R["comptage"]["prises_type"]) / (
    R["comptage"]["depots"] + R["comptage"]["prises"]) - 1
print("comptage :", R["comptage"])

S["gdl_t"] = np.array([(t - JOURS[0]).total_seconds() / 3600 for t in inst])
S["gdl_x"] = velos; S["gdl_libres"] = libres
S["nat_t"] = np.array([(t - JOURS[0]).total_seconds() / 3600 for t in inst_n])
S["nat_x"] = velos_n
i_w, v_w, l_w = WE["gare_de_lyon"]
S["we_t"] = np.array([(t - datetime(2026, 5, 24)).total_seconds() / 3600 for t in i_w])
S["we_x"] = v_w

# =====================================================================
# 2. Estimations sur quelques fenêtres
# =====================================================================
R["fenetres"] = {}


def fenetre(cle, ii, vv, ll, a, b):
    idx = D.indices(ii, a, b)
    n_d, t_d, n_p, t_p = V.estimer_taux(vv[idx], ll[idx])
    R["fenetres"][cle] = dict(n_d=int(n_d), t_d=float(t_d), n_p=int(n_p), t_p=float(t_p),
                              lam=ic(n_d, t_d), mu=ic(n_p, t_p), rho=ic_rho(n_d, t_d, n_p, t_p),
                              x0=int(vv[idx[0]]), x1=int(vv[idx[-1]]))
    print(cle, R["fenetres"][cle])
    return idx


fenetre("gdl_mar_6_8", inst, velos, libres, datetime(2026, 6, 2, 6), datetime(2026, 6, 2, 8))
fenetre("gdl_mer_6_8", inst, velos, libres, datetime(2026, 6, 3, 6), datetime(2026, 6, 3, 8))
fenetre("gdl_mar_9_19", inst, velos, libres, datetime(2026, 6, 2, 9), datetime(2026, 6, 2, 19))
fenetre("gdl_mar_18_20", inst, velos, libres, datetime(2026, 6, 2, 18), datetime(2026, 6, 2, 20))
fenetre("gdl_semaine", inst, velos, libres, inst[0], inst[-1])
fenetre("nat_mar_9_19", inst_n, velos_n, libres_n, datetime(2026, 6, 2, 9), datetime(2026, 6, 2, 19))
fenetre("nat_semaine", inst_n, velos_n, libres_n, inst_n[0], inst_n[-1])
fenetre("gdl_dimanche", i_w, v_w, l_w, i_w[0], i_w[-1])

# =====================================================================
# 3. Forçage moyen (Gare de Lyon, 4 jours ouvrés)
# =====================================================================
CPT = {j: comptes_creneaux(inst, velos, libres, j) for j in JOURS}
s_all = sum(CPT.values())
lam_all, mu_all = taux_depuis(s_all)
forc = []
for k in range(NB):
    l_, m_ = lam_all[k], mu_all[k]
    ligne = dict(n_d=int(s_all[k, 0]), t_d=float(s_all[k, 1]), n_p=int(s_all[k, 2]), t_p=float(s_all[k, 3]),
                 lam=ic(s_all[k, 0], s_all[k, 1]), mu=ic(s_all[k, 2], s_all[k, 3]))
    if l_ > 0 and m_ > 0:
        pi = V.loi_stationnaire(C, l_, m_)
        ligne.update(rho=float(l_ / m_), relax=float(1 / V.trou_spectral(C, l_, m_)),
                     melange=t_melange(C, l_, m_), moy_stat=float(pi @ np.arange(C + 1)))
    forc.append(ligne)
R["forcage"] = forc
R["forcage_derive_jour"] = float(((lam_all - mu_all) * HC).sum())
R["forcage_melange_min"] = min(f["melange"] for f in forc if "melange" in f)
R["forcage_melange_max"] = max(f["melange"] for f in forc if "melange" in f)
R["forcage_relax_min"] = min(f["relax"] for f in forc if "relax" in f)
R["forcage_relax_max"] = max(f["relax"] for f in forc if "relax" in f)
jour_ok = [f for f in forc[4:11] if "relax" in f]      # créneaux 8h-22h
R["forcage_relax_jour"] = [min(f["relax"] for f in jour_ok), max(f["relax"] for f in jour_ok)]
R["forcage_melange_jour"] = [min(f["melange"] for f in jour_ok), max(f["melange"] for f in jour_ok)]
S["lam_all"] = lam_all; S["mu_all"] = mu_all
S["lam_ic"] = np.array([[f["lam"][1], f["lam"][2]] for f in forc])
S["mu_ic"] = np.array([[f["mu"][1], f["mu"][2]] for f in forc])
for f in forc:
    print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in f.items()})

# tests du rapport de vraisemblance
R["rv"] = dict(
    depots_jour=rv(s_all[:, 0], s_all[:, 1]), prises_jour=rv(s_all[:, 2], s_all[:, 3]),
    depots_8_20=rv(s_all[4:10, 0], s_all[4:10, 1]), prises_8_20=rv(s_all[4:10, 2], s_all[4:10, 3]))
for a, b, nom in [(0, 1, "depots"), (2, 3, "prises")]:
    st_, dd_ = 0.0, 0
    for k in range(NB):
        N = np.array([CPT[j][k, a] for j in JOURS]); T = np.array([CPT[j][k, b] for j in JOURS])
        if N.sum() == 0:
            continue
        s_, d_, _ = rv(N, T)
        st_ += s_; dd_ += d_
    R["rv"][nom + "_entre_jours"] = (float(st_), int(dd_), float(chi2.sf(st_, dd_)))
print("RV :", R["rv"])

# dimanche (Gare de Lyon)
s_dim = comptes_creneaux(i_w, v_w, l_w, datetime(2026, 5, 24))
R["dimanche"] = dict(n_d=s_dim[:, 0].tolist(), t_d=s_dim[:, 1].tolist(), n_p=s_dim[:, 2].tolist(), t_p=s_dim[:, 3].tolist())

# =====================================================================
# 4. Loi stationnaire, relaxation
# =====================================================================
def stationnaire(cle, ii, vv, ll, cap, a, b):
    idx = D.indices(ii, a, b)
    v, l = vv[idx], ll[idx]
    n_d, t_d, n_p, t_p = V.estimer_taux(v, l)
    la, m_ = n_d / t_d, n_p / t_p
    emp = np.bincount(v, minlength=cap + 1) / len(v)
    pi = V.loi_stationnaire(cap, la, m_)
    heures = (len(v) - 1) / 60
    lois = V.lois_par_minute(cap, [la], [m_], heures, int(v[0]))
    moy, q05, q95 = bande(lois)
    R.setdefault("stationnaire", {})[cle] = dict(
        lam=float(la), mu=float(m_), rho=float(la / m_), vt=vt(emp, pi), moy_stat=float(pi @ np.arange(cap + 1)),
        moy_emp=float(v.mean()), relax=float(1 / V.trou_spectral(cap, la, m_)), melange=t_melange(cap, la, m_),
        couverture=float(((v >= q05) & (v <= q95)).mean()), x0=int(v[0]), x1=int(v[-1]),
        moy_fin=float(moy[-1]), q05_fin=int(q05[-1]), q95_fin=int(q95[-1]))
    print(cle, R["stationnaire"][cle])
    return v, emp, pi, moy, q05, q95


v, emp, pi, moy, q05, q95 = stationnaire("gdl", inst, velos, libres, C, datetime(2026, 6, 2, 9), datetime(2026, 6, 2, 19))
S["st_x"] = v; S["st_emp"] = emp; S["st_pi"] = pi; S["st_moy"] = moy; S["st_q05"] = q05; S["st_q95"] = q95
stationnaire("nation", inst_n, velos_n, libres_n, Cn, datetime(2026, 6, 2, 9), datetime(2026, 6, 2, 19))

# décroissance de la distance à l'équilibre pour deux créneaux
tt = np.linspace(0, 48, 481)
for cle, k, n0 in [("equil", 4, 0), ("pointe", 3, 0)]:
    Q = V.generateur(C, lam_all[k], mu_all[k])
    pi_k = V.loi_stationnaire(C, lam_all[k], mu_all[k])
    P1 = expm(0.1 * Q)
    p = np.zeros(C + 1); p[n0] = 1.0
    d = []
    for _ in tt:
        d.append(vt(p, pi_k))
        p = p @ P1
    S[f"dist_{cle}"] = np.array(d)
    gam = V.trou_spectral(C, lam_all[k], mu_all[k])
    S[f"borne_{cle}"] = 0.5 * np.sqrt((1 - pi_k[n0]) / pi_k[n0]) * np.exp(-gam * tt)
    R.setdefault("dist", {})[cle] = dict(lam=float(lam_all[k]), mu=float(mu_all[k]), gamma=float(gam),
                                         d2h=float(d[20]), d24h=float(d[240]))
S["dist_t"] = tt
print("distances :", R["dist"])

# =====================================================================
# 5. Prévision à 2 h (hors échantillon)
# =====================================================================
def prevision_2h(ii, vv, ll, cap, jours, cpt):
    dans = tot = 0
    res = []
    for j in jours:
        autres = [x for x in jours if x != j]
        lam_o, mu_o = taux_depuis(sum(cpt[x] for x in autres))
        for k in range(NB):
            idx = D.indices(ii, j + timedelta(hours=HC * k), j + timedelta(hours=HC * (k + 1)))
            if len(idx) < 121:
                continue
            reel = vv[idx]
            lois = V.lois_par_minute(cap, [lam_o[k]], [mu_o[k]], HC, int(reel[0]))
            moy, q05, q95 = bande(lois)
            ok = (reel >= q05) & (reel <= q95)
            dans += ok[1:].sum(); tot += len(ok) - 1
            res.append((j, k, int(reel[0]), int(reel[-1]), float(moy[-1]), int(q05[-1]), int(q95[-1]), float(ok[1:].mean())))
    fin = np.array([r[3] for r in res]); deb = np.array([r[2] for r in res]); prev = np.array([r[4] for r in res])
    hors = sum(1 for r in res if not (r[5] <= r[3] <= r[6]))
    return dict(couverture=float(dans / tot), creneaux=len(res), hors_fin=int(hors),
                p_binom=float(binomtest(int(hors), len(res), 0.10).pvalue),
                rmse_modele=float(np.sqrt(((fin - prev) ** 2).mean())), rmse_naif=float(np.sqrt(((fin - deb) ** 2).mean())),
                mae_modele=float(np.abs(fin - prev).mean()), mae_naif=float(np.abs(fin - deb).mean())), res


R["prev2h"] = {}
R["prev2h"]["gdl"], res_gdl = prevision_2h(inst, velos, libres, C, JOURS, CPT)
CPTN = {j: comptes_creneaux(inst_n, velos_n, libres_n, j) for j in JN}
R["prev2h"]["nation"], res_nat = prevision_2h(inst_n, velos_n, libres_n, Cn, JN, CPTN)
print("prévision 2 h :", R["prev2h"])
pires = sorted(res_gdl, key=lambda r: r[7])[:6]
print("  pires créneaux GdL :", [(f"{r[0]:%a} {2*r[1]}h", r[2], r[3], round(r[4], 1), r[5], r[6], round(r[7], 2)) for r in pires])


def detail_2h(cle, jour, k):
    """Bande et trajectoires simulées pour un créneau (figure)."""
    autres = [x for x in JOURS if x != jour]
    lam_o, mu_o = taux_depuis(sum(CPT[x] for x in autres))
    idx = D.indices(inst, jour + timedelta(hours=HC * k), jour + timedelta(hours=HC * (k + 1)))
    reel = velos[idx]
    lois = V.lois_par_minute(C, [lam_o[k]], [mu_o[k]], HC, int(reel[0]))
    moy, q05, q95 = bande(lois)
    grille = np.arange(121) / 60
    sims = []
    for _ in range(12):
        T, E = V.gillespie(C, lam_o[k], mu_o[k], HC, int(reel[0]), rng)
        sims.append(V.echantillonner(T, E, grille))
    S[f"c_{cle}_x"] = reel; S[f"c_{cle}_moy"] = moy; S[f"c_{cle}_q05"] = q05; S[f"c_{cle}_q95"] = q95
    S[f"c_{cle}_sims"] = np.array(sims)
    R.setdefault("court", {})[cle] = dict(lam=float(lam_o[k]), mu=float(mu_o[k]), x0=int(reel[0]), x1=int(reel[-1]),
                                          moy_fin=float(moy[-1]), q05_fin=int(q05[-1]), q95_fin=int(q95[-1]),
                                          couverture=float(((reel >= q05) & (reel <= q95))[1:].mean()))


detail_2h("matin", JOURS[0], 3)     # mardi 6-8 h
detail_2h("soir", JOURS[0], 10)     # mardi 20-22 h
print("court :", R["court"])

# =====================================================================
# 6. Prévision à 24 h
# =====================================================================
def prevision_24h(ii, vv, ll, cap, jours, cpt, cle):
    lignes = {}
    for j in jours:
        autres = [x for x in jours if x != j]
        s_o = sum(cpt[x] for x in autres)
        lam_o, mu_o = taux_depuis(s_o)
        idx = D.indices(ii, j, j + timedelta(hours=24) - timedelta(minutes=1))
        reel = vv[idx]
        n0 = int(reel[0])
        lam_c = np.full(NB, s_o[:, 0].sum() / s_o[:, 1].sum())
        mu_c = np.full(NB, s_o[:, 2].sum() / s_o[:, 3].sum())
        lam_i, mu_i = taux_depuis(cpt[j])
        d = {}
        for nom, (la, m_) in {"constant": (lam_c, mu_c), "creneaux": (lam_o, mu_o), "interne": (lam_i, mu_i)}.items():
            lois = V.lois_par_minute(cap, la, m_, HC, n0)[:len(reel)]
            moy, q05, q95 = bande(lois)
            d[nom] = dict(mae=float(np.abs(moy - reel).mean()), couverture=float(((reel >= q05) & (reel <= q95)).mean()))
            S[f"j_{cle}_{j:%d}_{nom}_moy"] = moy; S[f"j_{cle}_{j:%d}_{nom}_q05"] = q05; S[f"j_{cle}_{j:%d}_{nom}_q95"] = q95
        d["n0"] = n0
        d["lam_c"] = float(lam_c[0]); d["mu_c"] = float(mu_c[0])
        S[f"j_{cle}_{j:%d}_x"] = reel
        lignes[f"{j:%d}"] = d
    for nom in ("constant", "creneaux", "interne"):
        lignes[f"moy_{nom}"] = dict(mae=float(np.mean([lignes[f"{j:%d}"][nom]["mae"] for j in jours])),
                                   couverture=float(np.mean([lignes[f"{j:%d}"][nom]["couverture"] for j in jours])))
    return lignes


R["prev24h"] = dict(gdl=prevision_24h(inst, velos, libres, C, JOURS, CPT, "gdl"),
                    nation=prevision_24h(inst_n, velos_n, libres_n, Cn, JN, CPTN, "nat"))
for st in R["prev24h"]:
    for k_, d in R["prev24h"][st].items():
        print(st, k_, d)

# =====================================================================
# 7. Régime périodique
# =====================================================================
M = np.eye(C + 1)
for k in range(NB):
    M = M @ expm(HC * V.generateur(C, lam_all[k], mu_all[k]))
vp = np.linalg.eigvals(M)
modules = np.sort(np.abs(vp))[::-1]
w, vec = np.linalg.eig(M.T)
nu = np.real(vec[:, np.argmax(np.real(w))]); nu = nu / nu.sum()
assert np.abs(nu @ M - nu).max() < 1e-10 and nu.min() > -1e-12
etats = np.arange(C + 1)
# lois nu_t à chaque minute de la journée
lois_p = [nu]
p = nu.copy()
for k in range(NB):
    P1 = expm(V.generateur(C, lam_all[k], mu_all[k]) / 60)
    for _ in range(120):
        p = p @ P1
        lois_p.append(p)
lois_p = np.array(lois_p)
assert np.abs(lois_p[-1] - nu).max() < 1e-9
moy_p, q05_p, q95_p = bande(lois_p)
lam_min = np.repeat(lam_all, 120); mu_min = np.repeat(mu_all, 120)
plein = lois_p[:-1, C]; vide = lois_p[:-1, 0]
R["periodique"] = dict(
    vp2=float(modules[1]), vp3=float(modules[2]), doeblin=float(M.min(axis=0).sum()),
    moy_minuit=float(nu @ etats), p_inf10=float(nu[:11].sum()), p_sup40=float(nu[40:].sum()), p_inf1=float(nu[:2].sum()),
    q05_minuit=int(V.quantile(nu, 0.05)), q95_minuit=int(V.quantile(nu, 0.95)),
    temps_plein=float(plein.mean()), temps_vide=float(vide.mean()),
    depots_refuses=float((lam_min * plein).sum() / lam_min.sum()),
    prises_refusees=float((mu_min * vide).sum() / mu_min.sum()),
    moy_min=float(moy_p.min()), moy_max=float(moy_p.max()),
    h_moy_min=float(np.argmin(moy_p) / 60), h_moy_max=float(np.argmax(moy_p) / 60))
S["per_moy"] = moy_p; S["per_q05"] = q05_p; S["per_q95"] = q95_p; S["nu"] = nu
print("périodique :", R["periodique"])

# quatre jours depuis mardi 0 h
idx = D.indices(inst, JOURS[0], JOURS[3] + timedelta(hours=24))
reel4 = velos[idx]
lois4 = V.lois_par_minute(C, np.tile(lam_all, 4), np.tile(mu_all, 4), HC, int(reel4[0]))[:len(reel4)]
moy4, q05_4, q95_4 = bande(lois4)
R["quatre_jours"] = dict(n0=int(reel4[0]), couverture=float(((reel4 >= q05_4) & (reel4 <= q95_4)).mean()),
                         mae=float(np.abs(moy4 - reel4).mean()),
                         moy_24=float(moy4[1440]), moy_48=float(moy4[2880]), reel_24=int(reel4[1440]), reel_48=int(reel4[2880]),
                         vt_48=vt(lois4[2880], nu), vt_24=vt(lois4[1440], nu), vt_72=vt(lois4[4320], nu))
S["q4_x"] = reel4; S["q4_moy"] = moy4; S["q4_q05"] = q05_4; S["q4_q95"] = q95_4
print("quatre jours :", R["quatre_jours"])

# Monte-Carlo : contrôle croisé avec Kolmogorov et fonctionnelles de trajectoire
NSIM = 4000
grille = np.arange(1441) / 60
n0s = rng.choice(C + 1, size=NSIM, p=nu / nu.sum())
sims = np.zeros((NSIM, 1441), dtype=int)
for i in range(NSIM):
    T, E = V.gillespie_creneaux(C, lam_all, mu_all, HC, int(n0s[i]), rng)
    sims[i] = V.echantillonner(T, E, grille)
moy_mc = sims.mean(axis=0)
ecart_type = sims.std(axis=0).max()
duree_pleine = (sims[:, :-1] == C).sum(axis=1)          # minutes à saturation dans la journée
R["mc"] = dict(nsim=NSIM, ecart_max=float(np.abs(moy_mc - moy_p).max()), erreur_type=float(ecart_type / np.sqrt(NSIM)),
               ecart_type_max=float(ecart_type),
               p_sature=float((duree_pleine > 0).mean()), duree_moy=float(duree_pleine.mean()),
               duree_med=float(np.median(duree_pleine)), duree_q05=float(np.percentile(duree_pleine, 5)),
               duree_q95=float(np.percentile(duree_pleine, 95)), duree_q25=float(np.percentile(duree_pleine, 25)),
               duree_q75=float(np.percentile(duree_pleine, 75)))
obs = []
for j in JOURS:
    idx = D.indices(inst, j, j + timedelta(hours=24) - timedelta(minutes=1))
    obs.append(int((libres[idx] == 0).sum()))
R["mc"]["duree_obs"] = obs
S["mc_exemples"] = sims[:3]
print("MC :", R["mc"])

# =====================================================================
# 8. Variante « taux figés » et comptage unitaire (mardi, taux du jour)
# =====================================================================
def gillespie_fige(C, lam_tab, mu_tab, hc, duree, n0, rng):
    t, n = 0.0, n0
    T, E = [0.0], [n0]
    while t < duree:
        k = min(int(t // hc), len(lam_tab) - 1)
        lam = lam_tab[k] if n < C else 0.0
        mu = mu_tab[k] if n > 0 else 0.0
        tot = lam + mu
        if tot == 0:
            break
        t += rng.exponential(1 / tot)
        if t >= duree:
            break
        n += 1 if rng.random() < lam / tot else -1
        T.append(t); E.append(n)
    return np.array(T), np.array(E)


idx = D.indices(inst, JOURS[0], JOURS[0] + timedelta(hours=24) - timedelta(minutes=1))
reel = velos[idx]; lib = libres[idx]
n0 = int(reel[0])
lam_j, mu_j = taux_depuis(CPT[JOURS[0]])
# comptage unitaire (variations de +-1 seulement)
lam_u, mu_u = np.zeros(NB), np.zeros(NB)
for k in range(NB):
    a, b = 120 * k, min(120 * (k + 1), len(reel) - 1)
    d = np.diff(reel[a:b + 1])
    lam_u[k] = (d == 1).sum() / ((lib[a:b] > 0).sum() / 60)
    mu_u[k] = (d == -1).sum() / ((reel[a:b] > 0).sum() / 60)
grille = np.arange(len(reel)) / 60


def mae_exact(la, m_):
    lois = V.lois_par_minute(C, la, m_, HC, n0)[:len(reel)]
    moy, _, _ = bande(lois)
    return float(np.abs(moy - reel).mean()), float(moy[720])


def mae_fige(la, m_, nsim=3000):
    la, m_ = np.maximum(la, 0.3), np.maximum(m_, 0.3)
    tot = np.zeros(len(reel))
    for _ in range(nsim):
        T, E = gillespie_fige(C, la, m_, HC, 24.0, n0, rng)
        tot += V.echantillonner(T, E, grille)
    moy = tot / nsim
    return float(np.abs(moy - reel).mean()), float(moy[720])


R["variantes"] = dict(exact=mae_exact(lam_j, mu_j), fige=mae_fige(lam_j, mu_j),
                      unit_exact=mae_exact(lam_u, mu_u), unit_fige=mae_fige(lam_u, mu_u),
                      reel_12h=int(reel[720]), taux_nuit=float(lam_j[:3].sum() + mu_j[:3].sum()))
print("variantes :", R["variantes"])

# =====================================================================
# 9. Écarts à l'hypothèse de Poisson : arrivées groupées, rafales
# =====================================================================
def lots(ii, vv, ll, jours):
    obs_p = obs_m = 0
    att_p = att_m = 0.0
    for j in jours:
        for k in range(NB):
            idx = D.indices(ii, j + timedelta(hours=HC * k), j + timedelta(hours=HC * (k + 1)))
            if len(idx) < 2:
                continue
            n_d, t_d, n_p, t_p = V.estimer_taux(vv[idx], ll[idx])
            d = np.diff(vv[idx])
            obs_p += int((d >= 2).sum()); obs_m += int((d <= -2).sum())
            a = n_d / t_d / 60 if t_d > 0 else 0.0
            b = n_p / t_p / 60 if t_p > 0 else 0.0
            if a > 0 and b > 0:
                pp, pm = skellam.sf(1, a, b), skellam.cdf(-2, a, b)
            elif a > 0:
                pp, pm = poisson.sf(1, a), 0.0
            elif b > 0:
                pp, pm = 0.0, poisson.sf(1, b)
            else:
                pp = pm = 0.0
            att_p += pp * len(d); att_m += pm * len(d)
    return [obs_p, float(att_p), obs_m, float(att_m)]


dim = [datetime(2026, 5, 23), datetime(2026, 5, 24)]
R["lots"] = dict(
    gdl_semaine=lots(inst, velos, libres, JOURS), nation_semaine=lots(inst_n, velos_n, libres_n, JN),
    gdl_we=lots(*WE["gare_de_lyon"], dim), nation_we=lots(*WE["nation"], dim),
    princes_we=lots(*WE["parc_des_princes"], dim), terrasse_we=lots(*WE["claude_terrasse"], dim))
tot = np.array(list(R["lots"].values())).sum(axis=0)
R["lots"]["total"] = tot.tolist()
R["lots"]["rapport"] = float((tot[0] + tot[2]) / (tot[1] + tot[3]))
R["lots"]["p_total"] = float(poisson.sf(tot[0] + tot[2] - 1, tot[1] + tot[3]))
print("lots :", R["lots"])

# rafales
d_n = np.diff(velos_n)
i0 = inst_n.index(datetime(2026, 6, 3, 21, 36)); i1 = inst_n.index(datetime(2026, 6, 3, 21, 46))
s_mer = CPTN[JN[1]][10]; s_mar = CPTN[JN[0]][10]
R["rafales"] = dict(
    nation=dict(x0=int(velos_n[i0]), x1=int(velos_n[i1]), hausse=int(velos_n[i1] - velos_n[i0]),
                taux_creneau=float(s_mer[0] / s_mer[1]), taux_veille=float(s_mar[0] / s_mar[1]),
                p_creneau=float(poisson.sf(velos_n[i1] - velos_n[i0] - 1, s_mer[0] / s_mer[1] / 6)),
                p_veille=float(poisson.sf(velos_n[i1] - velos_n[i0] - 1, s_mar[0] / s_mar[1] / 6))))
j0 = i_w.index(datetime(2026, 5, 24, 20, 45)); j1 = i_w.index(datetime(2026, 5, 24, 20, 55))
jm = i_w.index(datetime(2026, 5, 24, 20, 53))
idx = D.indices(i_w, datetime(2026, 5, 24, 18), datetime(2026, 5, 24, 20))
n_d, t_d, n_p, t_p = V.estimer_taux(v_w[idx], l_w[idx])
R["rafales"]["gdl"] = dict(x0=int(v_w[j0]), x1=int(v_w[j1]), baisse=int(v_w[j0] - v_w[j1]),
                           saut_minute=int(v_w[jm] - v_w[jm + 1]), taux_18_20=float(n_p / t_p),
                           p_minute=float(poisson.sf(v_w[jm] - v_w[jm + 1] - 1, n_p / t_p / 60)),
                           p_10min=float(poisson.sf(v_w[j0] - v_w[j1] - 1, n_p / t_p / 6)))
print("rafales :", R["rafales"])

# =====================================================================
# 10. Exemples numériques pour la partie théorique
# =====================================================================
ex = {}
Q = V.generateur(C, lam_all[3], mu_all[3])                 # pointe du matin, 6-8 h
ex["pointe"] = dict(lam=float(lam_all[3]), mu=float(mu_all[3]),
                    sat_depuis_10=float(V.temps_moyen_atteinte(Q, [C])[10]),
                    sat_depuis_0=float(V.temps_moyen_atteinte(Q, [C])[0]))
Q = V.generateur(C, lam_all[4], mu_all[4])                 # creneau equilibre, 8-10 h
rho = lam_all[4] / mu_all[4]
ex["equil"] = dict(lam=float(lam_all[4]), mu=float(mu_all[4]), rho=float(rho),
                   sortie_depuis_30=float(V.temps_moyen_atteinte(Q, [0, C])[30]),
                   vide_depuis_30=float(V.temps_moyen_atteinte(Q, [0])[30]),
                   sat_depuis_30=float(V.temps_moyen_atteinte(Q, [C])[30]),
                   p_sat_avant_vide_30=float((1 - rho ** (-30)) / (1 - rho ** (-C))))
ex["relax_rho1"] = float(1 / V.trou_spectral(45, 4.0, 4.0))
ex["relax_rho1_approx"] = float(46 ** 2 / (np.pi ** 2 * 4.0))
ex["melange_rho1"] = t_melange(45, 4.0, 4.0)
R["exemples"] = ex
print("exemples :", ex)

# paires depot + prise dans la meme minute, invisibles dans le total (estimation)
paires = 0.0
for j in JOURS:
    l_, m_ = taux_depuis(CPT[j])
    paires += float((l_ * m_).sum() * HC / 60)
R["comptage"]["paires_attendues"] = paires
print("paires attendues :", paires)

with open(f"{OUT}/resultats.json", "w") as f:
    json.dump(R, f, indent=1, ensure_ascii=False)
np.savez(f"{OUT}/series.npz", **S)
print("OK")
