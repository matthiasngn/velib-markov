"""Assemble l'article : insère les valeurs calculées et le code dans les sources LaTeX.

Usage : python3 -I remplir.py <resultats.json> <dossier_src> <velib.py> <article.tex>
"""
import json
import math
import re
import sys

R = json.load(open(sys.argv[1]))
SRC, CODE, SORTIE = sys.argv[2], sys.argv[3], sys.argv[4]


def fr(x, nd=1):
    return f"{x:.{nd}f}".replace(".", "{,}")


def ent(x):
    return str(int(round(x)))


def adapt(x):
    if x >= 100:
        return ent(x)
    if x >= 10:
        return fr(x, 1)
    return fr(x, 2)


def pval(p):
    if p >= 0.01:
        return fr(p, 3 if p < 0.1 else 2)
    e = math.floor(math.log10(p))
    m = p / 10 ** e
    if p < 1e-6:
        m = round(m)
        if m == 10:
            m, e = 1, e + 1
        return f"{m}\\cdot10^{{{e}}}"
    return f"{fr(m, 1)}\\cdot10^{{{e}}}"


K = {}
fo = R["forcage"]
pe = R["periodique"]
va = R["variantes"]
st = R["stationnaire"]
fe = R["fenetres"]
p2 = R["prev2h"]
p24 = R["prev24h"]
mc = R["mc"]
ex = R["exemples"]
co = R["comptage"]
lo = R["lots"]
ra = R["rafales"]

# résumé
pmax = max(R["rv"]["depots_jour"][2], R["rv"]["prises_jour"][2])
K["pexp"] = str(math.floor(math.log10(pmax)) + 1)
K["couv2h_en"] = f"{100 * p2['gdl']['couverture']:.0f}"
K["melmin_ent"] = ent(R["forcage_melange_min"])
K["melmax_ent"] = ent(R["forcage_melange_max"])
K["vp2_en"] = f"{pe['vp2']:.2f}"
K["biais12"] = ent(va["exact"][1] - va["unit_fige"][1])

# section 3
K["relax_rho1"] = fr(ex["relax_rho1"], 1)
K["lam_pointe"] = fr(ex["pointe"]["lam"], 2)
K["mu_pointe"] = fr(ex["pointe"]["mu"], 2)
K["sat10"] = fr(ex["pointe"]["sat_depuis_10"], 1)
K["rho_equil"] = fr(ex["equil"]["rho"], 2)
K["sortie30"] = ent(ex["equil"]["sortie_depuis_30"])
K["psat30"] = fr(ex["equil"]["p_sat_avant_vide_30"], 2)

# sections 4 et 5
K["paires"] = ent(co["paires_attendues"])
assert co["depots_type"] - co["depots"] == co["prises_type"] - co["prises"]
K["paires_vues"] = str(co["depots_type"] - co["depots"])
K["nsim"] = f"{mc['nsim'] // 1000}\\,{mc['nsim'] % 1000:03d}"
K["mc_ecart"] = fr(mc["ecart_max"], 2)
K["mc_et"] = fr(mc["erreur_type"], 2)
K["perte_dep"] = ent(100 * co["perte_unit_depots"])
K["perte_pri"] = ent(100 * co["perte_unit_prises"])
baisses = sorted([va["exact"][1] - va["fige"][1], va["exact"][1] - va["unit_exact"][1]])
K["baisse_min"], K["baisse_max"] = ent(baisses[0]), ent(baisses[1])
for cle, nom in [("exact", "exact"), ("fige", "fige"), ("unit_exact", "unit"), ("unit_fige", "uf")]:
    K[f"v_{nom}_mae"] = fr(va[cle][0], 1)
    K[f"v_{nom}_12"] = fr(va[cle][1], 1)
K["reel12"] = str(va["reel_12h"])

# section 6.1
rows = ["\\multicolumn{7}{@{}l}{\\emph{Jours ouvrés, 1\\textsuperscript{er}--5 juin}}\\\\"]
for cle, nom in [("gdl_semaine", "Gare de Lyon"), ("nation_semaine", "Nation"),
                 (None, None), ("gdl_we", "Gare de Lyon"), ("nation_we", "Nation"),
                 ("princes_we", "Parc des Princes"), ("terrasse_we", "Claude Terrasse")]:
    if cle is None:
        rows.append("\\multicolumn{7}{@{}l}{\\emph{Week-end, 23--24 mai}}\\\\")
        continue
    d = R["donnees"][cle]
    rows.append(f"\\quad {nom} & {d['capacite']} & {fr(d['heures'], 1)} & {d['depots']} & {d['prises']}"
                f" & {d['min_pleine']} & {d['min_vide']}\\\\")
K["TAB_DONNEES"] = "\n".join(rows)

# section 6.2
f9 = fe["gdl_mar_9_19"]
for nom, (v, a, b) in [("lam919", f9["lam"]), ("mu919", f9["mu"]), ("rho919", f9["rho"])]:
    K[nom], K[nom + "a"], K[nom + "b"] = fr(v, 2), fr(a, 2), fr(b, 2)
K["estat919"] = fr(st["gdl"]["moy_stat"], 1)
K["emp919"] = fr(st["gdl"]["moy_emp"], 1)
K["vt919"] = fr(st["gdl"]["vt"], 3)
K["x0919"] = str(st["gdl"]["x0"])
K["couv919"] = ent(100 * st["gdl"]["couverture"])
K["relax919"] = fr(st["gdl"]["relax"], 1)
K["mel919"] = fr(st["gdl"]["melange"], 1)
K["relaxnat"] = fr(st["nation"]["relax"], 1)
K["melnat"] = fr(st["nation"]["melange"], 1)

rows = ["\\multicolumn{6}{@{}l}{\\emph{Gare de Lyon}}\\\\"]
for cle, nom in [("gdl_mar_6_8", "mardi 6--8~h"), ("gdl_mer_6_8", "mercredi 6--8~h"),
                 ("gdl_mar_9_19", "mardi 9--19~h"), ("gdl_mar_18_20", "mardi 18--20~h"),
                 ("gdl_semaine", "1\\textsuperscript{er}--5 juin"), ("gdl_dimanche", "23--24 mai"),
                 (None, None), ("nat_mar_9_19", "mardi 9--19~h"), ("nat_semaine", "1\\textsuperscript{er}--3 juin")]:
    if cle is None:
        rows.append("\\multicolumn{6}{@{}l}{\\emph{Nation}}\\\\")
        continue
    f = fe[cle]
    r, a, b = f["rho"]
    rows.append(f"\\quad {nom} & {f['n_d']} & {f['n_p']} & {adapt(f['lam'][0])} & {adapt(f['mu'][0])}"
                f" & {adapt(r)} \\ [{adapt(a)}\\,;\\,{adapt(b)}]\\\\")
K["TAB_FENETRES"] = "\n".join(rows)

# section 6.3
K["nd_pointe"] = str(fo[3]["n_d"])
K["np_pointe"] = str(fo[3]["n_p"])
K["rho_pointe"] = ent(fo[3]["rho"])
dim = R["dimanche"]
K["lam_dim68"] = fr(dim["n_d"][3] / dim["t_d"][3], 1)
K["lam_dim1012"] = fr(dim["n_d"][5] / dim["t_d"][5], 1)
rv = R["rv"]
K["w_dj"], K["p_dj"] = fr(rv["depots_jour"][0], 1), pval(rv["depots_jour"][2])
K["w_pj"], K["p_pj"] = fr(rv["prises_jour"][0], 1), pval(rv["prises_jour"][2])
K["p_d820"], K["p_p820"] = pval(rv["depots_8_20"][2]), pval(rv["prises_8_20"][2])
K["p_dej"], K["p_pej"] = pval(rv["depots_entre_jours"][2]), pval(rv["prises_entre_jours"][2])
assert rv["depots_jour"][1] == 11
K["melmin"] = fr(R["forcage_melange_min"], 1)
K["meljour_a"] = fr(R["forcage_melange_jour"][0], 1)
K["meljour_b"] = fr(R["forcage_melange_jour"][1], 1)
K["relax_pointe_min"] = ent(60 * fo[3]["relax"])
rows = []
for k, f in enumerate(fo):
    lab = f"{2 * k}--{2 * k + 2}~h"
    if "rho" in f:
        tr = fr(f["relax"], 2) if f["relax"] < 1 else fr(f["relax"], 1)
        rows.append(f"{lab} & {fr(f['lam'][0], 2)} & {fr(f['mu'][0], 2)} & {adapt(f['rho'])} & {tr} & {fr(f['melange'], 1)}\\\\")
    else:
        rows.append(f"{lab} & {fr(f['lam'][0], 2)} & {fr(f['mu'][0], 2)} & 0 & -- & --\\\\")
K["TAB_CRENEAUX"] = "\n".join(rows)

# section 6.4
g, n = p2["gdl"], p2["nation"]
K["n2h"] = str(g["creneaux"])
K["couv2h"] = fr(100 * g["couverture"], 1)
K["hors2h"] = str(g["hors_fin"])
K["att2h"] = fr(0.1 * g["creneaux"], 1)
K["pbin2h"] = pval(g["p_binom"])
K["rmse2h"], K["rmsenaif2h"] = fr(g["rmse_modele"], 1), fr(g["rmse_naif"], 1)
K["couv2hnat"] = fr(100 * n["couverture"], 1)
K["rmse2hnat"], K["rmsenaif2hnat"] = fr(n["rmse_modele"], 1), fr(n["rmse_naif"], 1)
so = R["court"]["soir"]
K["soir_x0"], K["soir_x1"], K["soir_prev"] = str(so["x0"]), str(so["x1"]), ent(so["moy_fin"])
gj = p24["gdl"]
K["eam24c"], K["eam24k"] = fr(gj["moy_constant"]["mae"], 1), fr(gj["moy_creneaux"]["mae"], 1)
K["couv24c"], K["couv24k"] = ent(100 * gj["moy_constant"]["couverture"]), ent(100 * gj["moy_creneaux"]["couverture"])
K["eam24i"] = fr(gj["moy_interne"]["mae"], 1)


def cellule(d):
    return f"{fr(d['mae'], 1)} ({ent(100 * d['couverture'])})"


rows = []
for st_, titre, jours in [("gdl", "Gare de Lyon", [("02", "mardi 2"), ("03", "mercredi 3"), ("04", "jeudi 4"), ("05", "vendredi 5")]),
                          ("nation", "Nation", [("02", "mardi 2"), ("03", "mercredi 3")])]:
    rows.append(f"\\multicolumn{{5}}{{@{{}}l}}{{\\emph{{{titre}}}}}\\\\")
    for j, nom in jours:
        d = p24[st_][j]
        rows.append(f"\\quad {nom} & {d['n0']} & {cellule(d['constant'])} & {cellule(d['creneaux'])} & {cellule(d['interne'])}\\\\")
    m = p24[st_]
    rows.append(f"\\quad moyenne & & {cellule(m['moy_constant'])} & {cellule(m['moy_creneaux'])} & {cellule(m['moy_interne'])}\\\\")
K["TAB_24H"] = "\n".join(rows)

# section 6.5
K["vp2"] = fr(pe["vp2"], 2)
K["doeblin"] = fr(pe["doeblin"], 2)
K["inv_vp2"] = ent(1 / pe["vp2"])
K["moyminuit"] = fr(pe["moy_minuit"], 1)
K["q05minuit"], K["q95minuit"] = str(pe["q05_minuit"]), str(pe["q95_minuit"])
K["pinf10"] = fr(pe["p_inf10"], 3)
q4 = R["quatre_jours"]
K["vt24"], K["vt48"], K["vt72"] = fr(q4["vt_24"], 2), fr(q4["vt_48"], 3), fr(q4["vt_72"], 3)
K["couv4j"] = fr(100 * q4["couverture"], 1)
K["tplein"] = fr(100 * pe["temps_plein"], 1)
K["deprefus"] = fr(100 * pe["depots_refuses"], 1)
K["tvide"] = fr(100 * pe["temps_vide"], 2)
K["psat"] = fr(100 * mc["p_sature"], 1)
K["dmed"], K["dq05"], K["dq95"] = ent(mc["duree_med"]), ent(mc["duree_q05"]), ent(mc["duree_q95"])
obs = [str(x) for x in mc["duree_obs"]]
K["dobs"] = ", ".join(obs[:-1]) + " et " + obs[-1]

# section 6.6
tot = lo["total"]
K["lots_obs"] = ent(tot[0] + tot[2])
K["lots_att"] = fr(tot[1] + tot[3], 1)
K["lots_rapport"] = fr(lo["rapport"], 1)
K["lots_p"] = pval(lo["p_total"])
rows = ["\\multicolumn{5}{@{}l}{\\emph{Jours ouvrés}}\\\\"]
for cle, nom in [("gdl_semaine", "Gare de Lyon"), ("nation_semaine", "Nation"), (None, None),
                 ("gdl_we", "Gare de Lyon"), ("nation_we", "Nation"), ("princes_we", "Parc des Princes"),
                 ("terrasse_we", "Claude Terrasse")]:
    if cle is None:
        rows.append("\\multicolumn{5}{@{}l}{\\emph{Week-end}}\\\\")
        continue
    a, b, c, d = lo[cle]
    rows.append(f"\\quad {nom} & {a} & {fr(b, 1)} & {c} & {fr(d, 1)}\\\\")
rows.append("\\midrule")
rows.append(f"Total & {ent(tot[0])} & {fr(tot[1], 1)} & {ent(tot[2])} & {fr(tot[3], 1)}\\\\")
K["TAB_LOTS"] = "\n".join(rows)
rn, rg = ra["nation"], ra["gdl"]
K["nat_x0"], K["nat_x1"], K["nat_hausse"] = str(rn["x0"]), str(rn["x1"]), str(rn["hausse"])
K["nat_taux_veille"] = fr(rn["taux_veille"], 1)
K["nat_p_veille"] = pval(rn["p_veille"])
K["gdl_saut"] = str(rg["saut_minute"])
K["gdl_p_minute"] = pval(rg["p_minute"])

# code
code = open(CODE).read().rstrip("\n")
assert all(ord(c) < 128 for c in code), "caractère non ASCII dans le code"
K["CODE"] = code

# assemblage
texte = "".join(open(f"{SRC}/partie{i}.tex").read() for i in (1, 2, 3))
manquants = sorted(set(re.findall(r"@@([A-Za-z0-9_]+)@@", texte)) - set(K))
assert not manquants, f"clés manquantes : {manquants}"
inutiles = sorted(set(K) - set(re.findall(r"@@([A-Za-z0-9_]+)@@", texte)))
texte = re.sub(r"@@([A-Za-z0-9_]+)@@", lambda m: K[m.group(1)], texte)
open(SORTIE, "w").write(texte)
print("clés inutilisées :", inutiles)
print("article écrit :", SORTIE, len(texte), "caractères")
