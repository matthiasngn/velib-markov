import numpy as np
from scipy.linalg import expm


# ---------- 1. Simulation ----------

def gillespie(C, lam, mu, duree, n0, rng):
    """Trajectoire a taux constants sur [0, duree]."""
    t = 0.0
    n = n0
    temps = [0.0]
    etats = [n0]
    while True:
        depot = lam if n < C else 0.0
        prise = mu if n > 0 else 0.0
        total = depot + prise
        if total == 0:
            break
        # delai avant le prochain evenement
        t = t + rng.exponential(1 / total)
        if t >= duree:
            break
        # nature de l'evenement
        if rng.random() < depot / total:
            n = n + 1
        else:
            n = n - 1
        temps.append(t)
        etats.append(n)
    return np.array(temps), np.array(etats)


def gillespie_creneaux(C, lam_tab, mu_tab, hc, n0, rng):
    """Taux constants sur des creneaux de duree hc."""
    temps = [0.0]
    etats = [n0]
    n = n0
    for k in range(len(lam_tab)):
        # simulation homogene sur le creneau k,
        # depuis l'etat atteint a la fin du precedent
        T, E = gillespie(C, lam_tab[k], mu_tab[k],
                         hc, n, rng)
        for j in range(1, len(T)):
            temps.append(k * hc + T[j])
            etats.append(E[j])
        n = E[-1]
    return np.array(temps), np.array(etats)


def echantillonner(temps, etats, grille):
    """Trajectoire lue aux instants de la grille."""
    valeurs = np.zeros(len(grille), dtype=int)
    j = 0
    for i in range(len(grille)):
        while (j + 1 < len(temps)
               and temps[j + 1] <= grille[i]):
            j = j + 1
        valeurs[i] = etats[j]
    return valeurs


# ---------- 2. Estimation ----------

def estimer_taux(velos, libres, pas=1 / 60):
    """Comptages et expositions (pas en heures)."""
    depots = 0
    prises = 0
    t_depot = 0.0   # temps ou un depot est possible
    t_prise = 0.0   # temps ou une prise est possible
    for i in range(len(velos) - 1):
        d = velos[i + 1] - velos[i]
        if d > 0:
            depots = depots + d
        if d < 0:
            prises = prises - d
        if libres[i] > 0:
            t_depot = t_depot + pas
        if velos[i] > 0:
            t_prise = t_prise + pas
    return depots, t_depot, prises, t_prise


def intervalle(nb, expo, z=1.96):
    """Taux estime et intervalle de confiance a 95 %."""
    taux = nb / expo
    marge = np.exp(z / np.sqrt(nb))
    return taux, taux / marge, taux * marge


# ---------- 3. Modele theorique ----------

def loi_stationnaire(C, lam, mu):
    rho = lam / mu
    poids = np.array([rho**n for n in range(C + 1)])
    return poids / poids.sum()


def generateur(C, lam, mu):
    Q = np.zeros((C + 1, C + 1))
    for n in range(C + 1):
        if n < C:
            Q[n, n + 1] = lam
        if n > 0:
            Q[n, n - 1] = mu
        Q[n, n] = -Q[n].sum()
    return Q


def trou_spectral(C, lam, mu):
    r = np.sqrt(lam * mu)
    return lam + mu - 2 * r * np.cos(np.pi / (C + 1))


def lois_par_minute(C, lam_tab, mu_tab, hc, n0):
    """Loi de X_t minute par minute (Kolmogorov)."""
    p = np.zeros(C + 1)
    p[n0] = 1.0
    lois = [p]
    for k in range(len(lam_tab)):
        Q = generateur(C, lam_tab[k], mu_tab[k])
        P = expm(Q / 60)      # transition sur une minute
        for i in range(round(60 * hc)):
            p = p @ P
            lois.append(p)
    return np.array(lois)


def quantile(loi, niveau):
    """Plus petit n tel que P(X <= n) >= niveau."""
    cumul = 0.0
    for n in range(len(loi)):
        cumul = cumul + loi[n]
        if cumul >= niveau:
            return n
    return len(loi) - 1


def temps_moyen_atteinte(Q, cible):
    """m(n) = E_n[T_cible] : Q m = -1 hors de la cible."""
    autres = [n for n in range(len(Q)) if n not in cible]
    A = Q[np.ix_(autres, autres)]
    m = np.zeros(len(Q))
    m[autres] = np.linalg.solve(A, -np.ones(len(autres)))
    return m
