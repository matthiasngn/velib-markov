"""Chargement des relevés sur une grille régulière d'une minute."""
import csv
from datetime import datetime, timedelta

import numpy as np

CAPACITE = {"gare_de_lyon": 45, "nation": 36,
            "parc_des_princes": 32, "claude_terrasse": 46}


def charger(dossier, fichiers):
    """Concatène des CSV et complète les minutes manquantes
    (report de la dernière valeur connue).
    Renvoie instants (liste de datetime), meca, elec, libres."""
    brut = {}
    for nom in fichiers:
        with open(f"{dossier}/{nom}") as f:
            for ligne in csv.DictReader(f):
                t = datetime.fromisoformat(ligne["instant"])
                brut[t] = (int(ligne["mecaniques"]),
                           int(ligne["electriques"]),
                           int(ligne["bornes_libres"]))
    t = min(brut)
    fin = max(brut)
    instants, valeurs = [], []
    courant = brut[t]
    while t <= fin:
        if t in brut:
            courant = brut[t]
        instants.append(t)
        valeurs.append(courant)
        t += timedelta(minutes=1)
    v = np.array(valeurs)
    return instants, v[:, 0], v[:, 1], v[:, 2]


def semaine_gdl(dossier):
    return charger(dossier, ["gare_de_lyon_20260601.csv",
                             "gare_de_lyon_20260603.csv"])


def semaine_nation(dossier):
    return charger(dossier, ["nation_20260601.csv"])


def weekend(dossier, station):
    return charger(dossier, [f"{station}_20260523.csv"])


def indices(instants, debut, fin):
    """Indices i tels que debut <= instants[i] <= fin."""
    return [i for i, t in enumerate(instants) if debut <= t <= fin]
