# Data

Minute-by-minute readings of four Vélib' stations in Paris, from the Vélib' Métropole real-time open data (© Smovengo / Syndicat Autolib' Vélib' Métropole), collected through [velib.nocle.fr](https://velib.nocle.fr).

## Files

Each file is one continuous recording, named after the station and the day it starts.

| File | First reading | Last reading | Rows |
|---|---|---|---|
| `gare_de_lyon_20260601.csv` | Mon 1 June 2026, 23:38 | Wed 3 June, 23:06 | 2,830 |
| `gare_de_lyon_20260603.csv` | Wed 3 June 2026, 23:29 | Fri 5 June, 23:06 | 2,838 |
| `nation_20260601.csv` | Mon 1 June 2026, 23:44 | Wed 3 June, 23:11 | 2,829 |
| `gare_de_lyon_20260523.csv` | Sat 23 May 2026, 20:56 | Sun 24 May, 20:55 | 1,432 |
| `nation_20260523.csv` | Sat 23 May 2026, 21:17 | Sun 24 May, 21:16 | 1,432 |
| `parc_des_princes_20260523.csv` | Sat 23 May 2026, 21:12 | Sun 24 May, 21:11 | 1,432 |
| `claude_terrasse_20260523.csv` | Sat 23 May 2026, 21:22 | Sun 24 May, 21:21 | 1,432 |

Columns:

| Column | Meaning |
|---|---|
| `instant` | time of the reading, local time (Europe/Paris), to the minute |
| `mecaniques` | mechanical bikes docked |
| `electriques` | electric bikes docked |
| `bornes_libres` | free docks |

The state of the model is $X_t$ = `mecaniques` + `electriques`. Capacities are 45 docks at Gare de Lyon, 36 at Nation, 32 at Parc des Princes and 46 at Claude Terrasse (`code/donnees.py`). One or two docks are often out of service, so $C - X_t$ can exceed the number of free docks.

## Cleaning

`code/parse_releves.py` turned the raw minute tables into these files: duplicate readings removed, recordings split wherever more than ten minutes are missing, and two misread capacity values discarded. The raw exports are not included.

About one minute per hour is missing. `code/donnees.py` fills each gap with the last known reading when it loads the files.
