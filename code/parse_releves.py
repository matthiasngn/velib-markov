"""Lecture des relevés Vélib (tableaux texte) -> CSV propres + diagnostic.

Usage : python3 -I parse_releves.py <dossier_source> <dossier_sortie>
"""
import csv
import os
import re
import sys
from collections import Counter
from datetime import datetime, timedelta

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

LIGNE = re.compile(
    r"\|\s*(\d{2}/\d{2}/\d{4} \d{2}:\d{2})\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*/\s*(\d+)\s*\|"
)

STATIONS = {
    "gare_de_lyon": [
        "Données_garde_de_lyon.docx",
        "Gare_de_lyon_2_juin.docx",
        "Gare_de_lyon_3_juin.docx",
        "Gare_de_lyon_4_juin.docx",
        "Gare_de_LYON_5_juin.docx",
    ],
    "nation": ["Données_nation.docx", "nation_2_juin.docx", "Nation_3_juin.docx"],
    "parc_des_princes": ["Données_par_des_princes.docx"],
    "claude_terrasse": ["Données_claude_terrasse.docx"],
}


def lire(chemin):
    lignes = []
    with open(chemin, encoding="utf-8") as f:
        for brut in f:
            m = LIGNE.search(brut)
            if m:
                t = datetime.strptime(m.group(1), "%d/%m/%Y %H:%M")
                lignes.append((t,) + tuple(int(m.group(k)) for k in (2, 3, 4, 5)))
    return lignes


for station, fichiers in STATIONS.items():
    print("=" * 70)
    print(station)
    releves = {}
    conflits = 0
    for nom in fichiers:
        rows = lire(os.path.join(SRC, nom))
        doublons_internes = len(rows) - len({r[0] for r in rows})
        for r in rows:
            if r[0] in releves and releves[r[0]] != r[1:]:
                conflits += 1
            releves[r[0]] = r[1:]
        print(f"  {nom:34s} {len(rows):5d} lignes  {rows[0][0]} -> {rows[-1][0]}"
              f"  doublons internes : {doublons_internes}")
    instants = sorted(releves)
    print(f"  total {len(instants)} relevés distincts, conflits entre fichiers : {conflits}")

    # champs "capacité" lus dans la colonne bornes libres x/C
    caps = Counter(releves[t][3] for t in instants)
    print("  valeurs du champ capacité :", dict(caps))
    cap = caps.most_common(1)[0][0]

    # découpage en segments continus (trou > 10 min => nouveau segment)
    segments, courant = [], [instants[0]]
    for a, b in zip(instants, instants[1:]):
        if b - a > timedelta(minutes=10):
            segments.append(courant)
            courant = []
        courant.append(b)
    segments.append(courant)

    for k, seg in enumerate(segments):
        manquants = int((seg[-1] - seg[0]).total_seconds() // 60) + 1 - len(seg)
        velos = [releves[t][0] + releves[t][1] for t in seg]
        libres = [releves[t][2] for t in seg]
        ecart = Counter(cap - v - l for v, l in zip(velos, libres))
        sauts = Counter(b - a for a, b in zip(velos, velos[1:]))
        print(f"  segment {k}: {seg[0]} ({seg[0]:%A}) -> {seg[-1]} ({seg[-1]:%A})  "
              f"{len(seg)} relevés, {manquants} minutes manquantes")
        print(f"     vélos min/max = {min(velos)}/{max(velos)} ; bornes libres min = {min(libres)}")
        print(f"     C - vélos - libres : {dict(sorted(ecart.items()))}")
        print(f"     variations par minute : {dict(sorted(sauts.items()))}")
        nom_csv = os.path.join(OUT, f"{station}_{seg[0]:%Y%m%d}.csv")
        with open(nom_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["instant", "mecaniques", "electriques", "bornes_libres"])
            for t in seg:
                w.writerow([t.isoformat(timespec="minutes"), *releves[t][:3]])
        print(f"     -> {nom_csv}  (capacité retenue C = {cap})")
