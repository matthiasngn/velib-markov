"""Numerical checks of velib.py against the closed-form results of the paper."""

import numpy as np
import pytest
from scipy.linalg import expm

import velib as V

CASES = [(10, 2.0, 3.0), (45, 10.05, 0.5), (36, 2.1, 2.5), (8, 4.0, 4.0)]


def tv(p, q):
    return 0.5 * np.abs(p - q).sum()


# ---------- stationary law (Proposition 3.1) ----------

@pytest.mark.parametrize("C, lam, mu", CASES)
def test_stationary_law_is_invariant(C, lam, mu):
    pi = V.loi_stationnaire(C, lam, mu)
    assert pi.sum() == pytest.approx(1.0)
    assert np.allclose(pi @ V.generateur(C, lam, mu), 0.0, atol=1e-12)


def test_stationary_law_closed_form():
    C, rho = 12, 0.7
    n = np.arange(C + 1)
    assert np.allclose(V.loi_stationnaire(C, rho, 1.0), rho**n * (1 - rho) / (1 - rho ** (C + 1)))
    assert np.allclose(V.loi_stationnaire(C, 3.0, 3.0), 1 / (C + 1))


# ---------- spectrum (Theorem 3.3) ----------

@pytest.mark.parametrize("C, lam, mu", CASES)
def test_spectrum_of_the_generator(C, lam, mu):
    k = np.arange(1, C + 1)
    s = -(lam + mu) + 2 * np.sqrt(lam * mu) * np.cos(k * np.pi / (C + 1))
    expected = np.sort(np.append(s, 0.0))
    computed = np.sort(np.linalg.eigvals(V.generateur(C, lam, mu)).real)
    assert np.allclose(computed, expected, atol=1e-8 * (lam + mu))


@pytest.mark.parametrize("C, lam, mu", CASES)
def test_eigenvectors(C, lam, mu):
    Q, rho, n = V.generateur(C, lam, mu), lam / mu, np.arange(C + 1)
    for k in (1, 2, C):
        theta = k * np.pi / (C + 1)
        f = rho ** (-n / 2) * (np.sin((n + 1) * theta) - rho ** (-0.5) * np.sin(n * theta))
        s = -(lam + mu) + 2 * np.sqrt(lam * mu) * np.cos(theta)
        assert np.allclose(Q @ f, s * f, atol=1e-8 * np.abs(f).max() * (lam + mu))


@pytest.mark.parametrize("C, lam, mu", CASES)
def test_spectral_gap(C, lam, mu):
    gap = V.trou_spectral(C, lam, mu)
    eig = np.sort(np.linalg.eigvals(V.generateur(C, lam, mu)).real)
    assert gap == pytest.approx(-eig[-2], rel=1e-8)
    closed = (np.sqrt(lam) - np.sqrt(mu)) ** 2 + 4 * np.sqrt(lam * mu) * np.sin(np.pi / (2 * (C + 1))) ** 2
    assert gap == pytest.approx(closed, rel=1e-10)


# ---------- convergence (Theorems 2.6 and 3.4) ----------

@pytest.mark.parametrize("C, lam, mu", [(10, 2.0, 3.0), (20, 4.0, 4.0)])
def test_distance_to_equilibrium_bounds(C, lam, mu):
    Q, pi, gap = V.generateur(C, lam, mu), V.loi_stationnaire(C, lam, mu), V.trou_spectral(C, lam, mu)
    rho, n, theta = lam / mu, np.arange(C + 1), np.pi / (C + 1)
    f1 = rho ** (-n / 2) * (np.sin((n + 1) * theta) - rho ** (-0.5) * np.sin(n * theta))
    for i in (0, C // 2, C):
        for t in (0.5, 2.0, 10.0):
            d = tv(expm(t * Q)[i], pi)
            lower = abs(f1[i]) / (2 * np.abs(f1).max()) * np.exp(-gap * t)
            upper = 0.5 * np.sqrt((1 - pi[i]) / pi[i]) * np.exp(-gap * t)
            assert lower - 1e-12 <= d <= upper + 1e-12


def test_doeblin_contraction():
    C, lam, mu, h = 10, 2.0, 3.0, 1.0
    Q, pi = V.generateur(C, lam, mu), V.loi_stationnaire(C, lam, mu)
    delta = expm(h * Q).min(axis=0).sum()
    assert delta > 0
    for t in (1.0, 3.0, 7.5):
        assert tv(expm(t * Q)[0], pi) <= (1 - delta) ** np.floor(t / h) + 1e-12


# ---------- hitting times (Corollary 3.6) ----------

def mean_time_to_empty(C, lam, mu):
    n, rho = np.arange(C + 1), lam / mu
    if rho == 1:
        return n * (2 * C - n + 1) / (2 * mu)
    return (n - rho ** (C + 1 - n) * (1 - rho**n) / (1 - rho)) / (mu * (1 - rho))


@pytest.mark.parametrize("lam, mu", [(2.0, 3.0), (4.0, 4.0), (6.0, 2.0)])
def test_mean_time_to_empty(lam, mu):
    C = 15
    m = V.temps_moyen_atteinte(V.generateur(C, lam, mu), [0])
    assert np.allclose(m, mean_time_to_empty(C, lam, mu), rtol=1e-8)


def test_mean_time_to_fill_at_the_morning_rush():
    # Section 3.4: Gare de Lyon, 6-8 am rates, a station holding 10 bikes fills up in 3.7 h on average
    C, lam, mu = 45, 10.050761421319809, 0.5
    m = V.temps_moyen_atteinte(V.generateur(C, lam, mu), [C])
    assert m[10] == pytest.approx(3.664629285144825, rel=1e-9)
    # filling is emptying for the mirrored station (bikes <-> free docks)
    assert np.allclose(m, mean_time_to_empty(C, mu, lam)[::-1], rtol=1e-8)


@pytest.mark.parametrize("lam, mu", [(2.0, 3.0), (4.0, 4.0)])
def test_probability_to_fill_before_emptying(lam, mu):
    C, rho = 12, lam / mu
    Q = V.generateur(C, lam, mu)
    inner = np.arange(1, C)
    h = np.zeros(C + 1)
    h[C] = 1.0
    h[inner] = np.linalg.solve(Q[np.ix_(inner, inner)], -Q[inner, C])  # Qh = 0 inside, h(0)=0, h(C)=1
    n = np.arange(C + 1)
    closed = n / C if rho == 1 else (1 - rho ** (-n)) / (1 - rho ** (-C))
    assert np.allclose(h, closed)


# ---------- exact law minute by minute ----------

def test_laws_per_minute_match_the_semigroup():
    C, n0 = 20, 5
    lam_tab, mu_tab = np.array([3.0, 0.5, 6.0]), np.array([1.0, 4.0, 6.0])
    laws = V.lois_par_minute(C, lam_tab, mu_tab, 2, n0)
    assert laws.shape == (3 * 120 + 1, C + 1)
    assert np.allclose(laws.sum(axis=1), 1.0)
    p = np.zeros(C + 1)
    p[n0] = 1.0
    for lam, mu in zip(lam_tab, mu_tab):
        p = p @ expm(2.0 * V.generateur(C, lam, mu))
    assert np.allclose(laws[-1], p, atol=1e-10)


def test_quantile():
    law = np.array([0.1, 0.2, 0.4, 0.3])
    assert V.quantile(law, 0.05) == 0
    assert V.quantile(law, 0.3) == 1
    assert V.quantile(law, 0.95) == 3


# ---------- simulation ----------

def occupation(times, states, C, t_end):
    out = np.zeros(C + 1)
    np.add.at(out, states, np.diff(np.append(times, t_end)))
    return out / out.sum()


def test_gillespie_settles_on_the_stationary_law():
    C, lam, mu, T = 15, 3.0, 3.5, 20_000.0
    times, states = V.gillespie(C, lam, mu, T, C // 2, np.random.default_rng(1))
    assert ((states >= 0) & (states <= C)).all()
    assert (np.abs(np.diff(states)) == 1).all()
    assert tv(occupation(times, states, C, T), V.loi_stationnaire(C, lam, mu)) < 0.02


def test_slots_with_equal_rates_behave_like_constant_rates():
    C, lam, mu = 10, 2.0, 2.5
    times, states = V.gillespie_creneaux(C, np.full(40, lam), np.full(40, mu), 500.0, 5, np.random.default_rng(2))
    assert tv(occupation(times, states, C, 40 * 500.0), V.loi_stationnaire(C, lam, mu)) < 0.02


def test_simulated_mean_matches_the_exact_law():
    C, n0, hc = 20, 2, 2.0
    lam_tab, mu_tab = np.array([6.0, 1.0, 3.0]), np.array([1.0, 5.0, 3.0])
    grid = np.arange(0, 3 * hc + 1e-9, 1 / 60)
    rng = np.random.default_rng(3)
    paths = np.array([
        V.echantillonner(*V.gillespie_creneaux(C, lam_tab, mu_tab, hc, n0, rng), grid) for _ in range(3000)
    ])
    exact = V.lois_par_minute(C, lam_tab, mu_tab, hc, n0) @ np.arange(C + 1)
    assert np.abs(paths.mean(axis=0) - exact).max() < 0.3


def test_echantillonner():
    times, states = np.array([0.0, 1.0, 2.5]), np.array([3, 4, 3])
    assert V.echantillonner(times, states, np.array([0.0, 0.99, 1.0, 2.0, 3.0])).tolist() == [3, 3, 4, 4, 3]


# ---------- estimation (Section 4) ----------

def test_counts_and_exposures():
    velos = np.array([2, 3, 3, 1, 0, 0, 2])
    libres = np.array([3, 2, 0, 4, 5, 5, 3])
    depots, t_depot, prises, t_prise = V.estimer_taux(velos, libres)
    assert (depots, prises) == (3, 3)  # +1, +2 / -2, -1
    assert t_depot == pytest.approx(5 / 60)  # minutes with a free dock before each step
    assert t_prise == pytest.approx(4 / 60)  # minutes with a bike before each step


def test_confidence_interval():
    rate, low, high = V.intervalle(100, 20.0)
    assert rate == pytest.approx(5.0)
    assert low == pytest.approx(5.0 * np.exp(-1.96 / 10))
    assert high == pytest.approx(5.0 * np.exp(1.96 / 10))


def test_rates_recovered_from_minute_readings():
    C, lam, mu = 60, 2.0, 2.0
    times, states = V.gillespie(C, lam, mu, 5_000.0, 30, np.random.default_rng(4))
    velos = V.echantillonner(times, states, np.arange(0.0, 5_000.0, 1 / 60))
    depots, t_depot, prises, t_prise = V.estimer_taux(velos, C - velos)
    # minute sampling hides a deposit and a pick-up falling in the same minute (Section 4.3)
    assert depots / t_depot == pytest.approx(lam, rel=0.08)
    assert prises / t_prise == pytest.approx(mu, rel=0.08)
    assert depots / t_depot - prises / t_prise == pytest.approx(lam - mu, abs=0.1)
