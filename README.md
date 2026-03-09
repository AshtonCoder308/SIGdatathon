# SIG Datathon — London Underground Network Analysis

Analysis of the TfL (Transport for London) Underground network, focusing on station-level passenger flows and crowding patterns across pre-COVID, recovery, and recent periods (2019, 2022, 2023).

## Project Overview

This project builds a graph representation of the London Underground and overlays it with real-world ridership data to explore how passenger behaviour has changed since COVID-19.

**Key questions:**
- How have station entry/exit counts recovered since the 2019 pre-COVID baseline?
- Which lines and stations are most/least congested?
- How does crowding vary by time of day and day of week?

## Data Sources

| Dataset | Source | Description |
|---------|--------|-------------|
| Station geometry | [Oliver O'Brien / OpenStreetMap](https://github.com/oobrien/vis) | GeoJSON with station coordinates and line branch polylines |
| Annual station counts | [TfL Crowding Data](https://crowding.data.tfl.gov.uk/Annual%20Station%20Counts) | Annualised entry/exit counts (2019, 2022, 2023) |
| Crowding averages | [TfL API](https://api-portal.tfl.gov.uk/) | Hourly busyness percentages per station per day of week |

## Setup

### Prerequisites

- Python 3.10+
- A [TfL API key](https://api-portal.tfl.gov.uk/) (free registration required for crowding data)

### Install dependencies

```bash
pip install requests pandas openpyxl tqdm
```

For graph construction (next steps):

```bash
pip install networkx
```

### Set your TfL API key

**Command Prompt:**
```cmd
set TFL_API_KEY=your_key_here
```

**PowerShell:**
```powershell
$env:TFL_API_KEY = "your_key_here"
```

Or set it as a permanent environment variable via Windows System Settings.

## Project Structure

```
SIG Datathon/
├── David/
│   ├── sig_datathon_david.ipynb   # Data ingestion & parsing notebook
│   └── data/                      # Downloaded datasets (gitignored)
│       ├── tfl_stations.json
│       ├── tfl_lines.json
│       ├── annual_counts_2019.xlsx
│       ├── annual_counts_2022.xlsx
│       ├── annual_counts_2023.xlsx
│       └── station_crowding.json
└── README.md
```

## Methodology

### 1. Data Ingestion (`David/sig_datathon_david.ipynb`)

**Station & Line Geometry**
Downloads GeoJSON from Oliver O'Brien's TubeCreature project, flattened into a DataFrame with one row per station (`stations_df`: 552 stations, columns: `station_id`, `name`, `zone`, `lines`, `lat`, `lon`).

**Graph Edge Construction**
Walks each line branch's coordinate polyline and snaps coordinates to the nearest station (within 150 m threshold). Station-to-station pairs become edges with:
- Haversine distance (km)
- Estimated travel time (minutes), assuming average tube speed of 33 km/h

Result: `edges_df` with 553 edges across 17 lines.

**Annual Entry/Exit Counts**
Parses TfL's Excel files (different column schemas per year) into a unified format with day-of-week entry/exit columns (`mon`, `twt`, `fri`, `sat`, optionally `sun`) and annualised totals. Covers:
- 2019 — pre-COVID peak baseline (425 stations)
- 2022 — COVID recovery year (427 stations)
- 2023 — most recent complete year (427 stations)

**Crowding Data**
Fetches per-station hourly busyness percentages from the TfL Crowding API across all 11 Underground lines (272 unique stations), stored in `station_crowding.json`.

### 2. Graph Construction (next step)

Build a NetworkX graph from `stations_df` and `edges_df`, annotated with ridership and crowding data, for network analysis and visualisation.

## Datasets Summary

| DataFrame | Shape | Description |
|-----------|-------|-------------|
| `stations_df` | 552 × 6 | One row per TfL station with coordinates and zone |
| `edges_df` | 553 × 7 | Station-to-station edges with distance and travel time |
| `annual_dfs['2019']` | 425 × 15 | Entry/exit counts, pre-COVID baseline |
| `annual_dfs['2022']` | 427 × 19 | Entry/exit counts, COVID recovery |
| `annual_dfs['2023']` | 427 × 19 | Entry/exit counts, most recent |
| `crowding_df` | 272 × 7 | Hourly busyness by station, day, and time band |

## Team

- **Ashton** — project setup, coordination
- **David** — data ingestion and parsing pipeline
