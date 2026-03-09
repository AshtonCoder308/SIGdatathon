"""
Spatial Demand Surface Model — Royal Borough of Greenwich
=========================================================
Datathon 2025/26 | LSESU Mathematics x Data Science

This script:
  1. Builds a population density field ρ(x,y) from Gaussian kernels
  2. Builds a transport accessibility field A(x,y) from 17 station nodes
  3. Computes the unmet demand surface D(x,y) = ρ(x,y) · (1 - A/Amax)
  4. Runs sensitivity analysis over (λ, weight scheme) — 20 combinations
  5. Generates 4 publication-quality figures

Requirements: numpy, scipy, matplotlib
"""

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from collections import Counter

# ============================================================
# 1. GRID SETUP
# ============================================================
# Borough bounds (approx): lat 51.43–51.50, lon −0.02 to 0.13
# Resolution: ~100m (0.001° lat ≈ 111m, 0.0015° lon ≈ 104m at 51.5°N)

LAT_MIN, LAT_MAX = 51.430, 51.500
LON_MIN, LON_MAX = -0.020, 0.130
LAT_STEP = 0.001
LON_STEP = 0.0015

lats = np.arange(LAT_MIN, LAT_MAX, LAT_STEP)
lons = np.arange(LON_MIN, LON_MAX, LON_STEP)
LON, LAT = np.meshgrid(lons, lats)

print(f"Grid: {len(lats)} × {len(lons)} = {len(lats)*len(lons)} points")
print(f"Lat range: [{lats[0]:.3f}, {lats[-1]:.3f}]")
print(f"Lon range: [{lons[0]:.4f}, {lons[-1]:.4f}]")


# ============================================================
# 2. STATION DATA
# ============================================================
# Each station: (latitude, longitude, base_weight, name, mode)
# Weights reflect capacity × frequency hierarchy (see report)
#   Tube: 8.0    | ~36,000 pax/hr (1200/train × 30 tph)
#   Elizabeth: 5.0| ~18,000 pax/hr (1500/train × 12 tph)
#   DLR: 4.0–5.5 | ~10,500 pax/hr (600/train × 15–20 tph)
#   Rail: 2.0–3.5| ~4,000 pax/hr  (800/train × 4–6 tph)

STATIONS = [
    # Tube
    (51.5005, 0.0039, 8.0, "North Greenwich", "tube"),
    # DLR
    (51.4827, -0.0106, 5.0, "Cutty Sark", "dlr"),
    (51.4791, -0.0155, 4.5, "Greenwich", "dlr"),
    (51.4882, 0.0698, 5.5, "Woolwich Arsenal", "dlr"),
    (51.4907, 0.0712, 5.0, "Woolwich", "dlr"),
    (51.4678, -0.0145, 3.0, "Lewisham", "dlr"),
    # Elizabeth Line
    (51.4910, 0.0799, 5.0, "Woolwich (EL)", "elizabeth"),
    (51.4906, 0.1204, 5.0, "Abbey Wood (EL)", "elizabeth"),
    # National Rail (Southeastern)
    (51.4905, 0.0087, 4.0, "Maze Hill", "rail"),
    (51.4950, 0.0380, 4.5, "Charlton", "rail"),
    (51.4625, 0.0282, 3.0, "Kidbrooke", "rail"),
    (51.4513, 0.0524, 3.5, "Eltham", "rail"),
    (51.4655, 0.0091, 3.0, "Blackheath", "rail"),
    (51.4563, 0.0714, 2.5, "Falconwood", "rail"),
    (51.4688, 0.0455, 2.5, "Lee", "rail"),
    (51.4720, 0.0648, 2.0, "Westcombe Pk", "rail"),
    (51.4903, 0.0955, 3.0, "Plumstead", "rail"),
]


# ============================================================
# 3. POPULATION DENSITY FIELD  ρ(x, y)
# ============================================================
# Modelled as superposition of Gaussian kernels centred on known
# residential/commercial centres, calibrated to ONS Census 2021.
#
# ρ(x,y) = ρ₀ + Σ wₖ · exp(−d²/(2σₖ²))
#
# Each centre: (lat, lon, weight, spread_metres)

POP_CENTRES = [
    (51.478, -0.010, 1.00, 600),   # Greenwich town centre
    (51.462,  0.026, 0.90, 700),   # Kidbrooke Village
    (51.451,  0.052, 0.85, 800),   # Eltham
    (51.490,  0.069, 0.80, 700),   # Woolwich
    (51.495,  0.038, 0.70, 600),   # Charlton
    (51.500,  0.004, 0.75, 500),   # Greenwich Peninsula
    (51.490,  0.120, 0.60, 800),   # Abbey Wood / Thamesmead
    (51.495,  0.100, 0.65, 900),   # Thamesmead
    (51.466,  0.009, 0.50, 500),   # Blackheath edge
    (51.457,  0.072, 0.60, 600),   # Falconwood / Welling edge
    (51.470,  0.050, 0.50, 500),   # Lee / Eltham north
    (51.485,  0.095, 0.55, 700),   # Plumstead
    (51.462,  0.040, 0.60, 500),   # Mottingham edge
]

# Negative kernels for parks / low-density areas
LOW_DENSITY = [
    (51.468, 0.002, -0.40, 400),   # Blackheath common
    (51.478, 0.000, -0.30, 300),   # Greenwich Park
    (51.458, 0.032, -0.20, 300),   # Sutcliffe Park
]

RHO_BASELINE = 0.3  # baseline density everywhere


def build_population_field():
    """Build ρ(x,y) on the grid."""
    rho = np.ones_like(LAT) * RHO_BASELINE

    for plat, plon, weight, spread in POP_CENTRES + LOW_DENSITY:
        dlat_m = (LAT - plat) * 111_000       # degrees → metres (lat)
        dlon_m = (LON - plon) * 69_000         # degrees → metres (lon at 51.5°N)
        dist_sq = dlat_m**2 + dlon_m**2
        rho += weight * np.exp(-dist_sq / (2 * spread**2))

    return np.clip(rho, 0.05, None)


# ============================================================
# 4. ACCESSIBILITY FIELD  A(x, y)
# ============================================================
# A(x,y) = Σⱼ wⱼ · exp(−d(x,y; xⱼ,yⱼ)² / 2λ²)
#
# λ = catchment decay parameter (metres)
# wⱼ = station weight (capacity × frequency index)


def build_accessibility_field(lam, weight_override=None):
    """
    Build A(x,y) on the grid.

    Parameters:
        lam : float — catchment decay in metres
        weight_override : dict or None — if provided, maps mode → weight
                          (overrides base_weight in STATIONS)
    Returns:
        accessibility : ndarray — raw accessibility values
        acc_norm : ndarray — normalised to [0, 1]
    """
    accessibility = np.zeros_like(LAT)

    for slat, slon, base_w, name, mode in STATIONS:
        if weight_override is not None:
            w = weight_override[mode]
        else:
            w = base_w

        dlat_m = (LAT - slat) * 111_000
        dlon_m = (LON - slon) * 69_000
        dist_sq = dlat_m**2 + dlon_m**2
        accessibility += w * np.exp(-dist_sq / (2 * lam**2))

    acc_norm = accessibility / accessibility.max()
    return accessibility, acc_norm


# ============================================================
# 5. DEMAND SURFACE  D(x, y)
# ============================================================
# D(x,y) = ρ(x,y) · (1 − A(x,y)/Amax)
# High where population is high AND accessibility is low.


def build_demand_surface(rho, acc_norm, smooth_sigma=2):
    """Compute unmet demand and apply Gaussian smoothing."""
    demand = rho * (1.0 - acc_norm)
    demand = gaussian_filter(demand, sigma=smooth_sigma)
    return demand


def find_peak(demand):
    """Return (lat, lon, value) of the global demand maximum."""
    idx = np.unravel_index(np.argmax(demand), demand.shape)
    return lats[idx[0]], lons[idx[1]], demand[idx]


# ============================================================
# 6. BUILD BASELINE MODEL
# ============================================================

print("\n=== Building baseline model (λ=800m, capacity-frequency weights) ===")

population = build_population_field()
_, acc_norm_baseline = build_accessibility_field(lam=800)
demand_baseline = build_demand_surface(population, acc_norm_baseline)

peak_lat, peak_lon, peak_val = find_peak(demand_baseline)
print(f"Peak unmet demand: ({peak_lat:.4f}°N, {peak_lon:.4f}°E)")
print(f"Peak demand value: {peak_val:.3f}")


# ============================================================
# 7. SENSITIVITY ANALYSIS
# ============================================================
# Vary: λ ∈ {400, 600, 800, 1000, 1200}
#       weight scheme ∈ {Equal, Mild, Baseline, Steep}
# → 20 combinations

LAMBDAS = [400, 600, 800, 1000, 1200]

WEIGHT_SCHEMES = {
    "Equal":    {"tube": 4.0, "elizabeth": 4.0, "dlr": 4.0, "rail": 4.0},
    "Mild":     {"tube": 6.0, "elizabeth": 5.0, "dlr": 4.0, "rail": 3.5},
    "Baseline": None,  # use base_weight from STATIONS
    "Steep":    {"tube": 12.0, "elizabeth": 8.0, "dlr": 4.0, "rail": 1.5},
}

SCHEME_LABELS = {
    "Equal":    "Equal (all 4.0)",
    "Mild":     "Mild (Tube 6, EL 5, DLR 4, Rail 3.5)",
    "Baseline": "Baseline (Tube 8, EL 5, DLR 4.5, Rail 2.5)",
    "Steep":    "Steep (Tube 12, EL 8, DLR 4, Rail 1.5)",
}


def classify_area(lat, lon):
    """Assign a named area label to peak coordinates."""
    if lat < 51.455 and lon > 0.04:
        return "Eltham"
    elif lat < 51.465 and lon > 0.02:
        return "Kidbrooke-Eltham"
    elif lat < 51.475 and lon < 0.02:
        return "Blackheath-Kidbrooke"
    elif lat < 51.485 and lon < 0.01:
        return "Greenwich TC"
    elif lat > 51.485 and lon > 0.08:
        return "Thamesmead"
    else:
        return "Kidbrooke corridor"


print("\n=== Sensitivity Analysis: 20 parameter combinations ===\n")
print(f"{'λ (m)':<8} {'Weight Scheme':<50} {'Peak Lat':>10} {'Peak Lon':>10} "
      f"{'Demand':>8} {'Winner':<20}")
print("=" * 110)

results = []
winner_grid = []  # for the figure

for lam in LAMBDAS:
    row = []
    for scheme_name in ["Equal", "Mild", "Baseline", "Steep"]:
        override = WEIGHT_SCHEMES[scheme_name]
        _, acc_n = build_accessibility_field(lam, weight_override=override)
        d = build_demand_surface(population, acc_n)
        plat, plon, pval = find_peak(d)
        area = classify_area(plat, plon)

        marker = " ★" if (lam == 800 and scheme_name == "Baseline") else ""
        print(f"{lam:<8} {SCHEME_LABELS[scheme_name]:<50} {plat:>10.4f} "
              f"{plon:>10.4f} {pval:>8.3f} {area:<20}{marker}")

        results.append({
            "lambda": lam,
            "scheme": scheme_name,
            "peak_lat": plat,
            "peak_lon": plon,
            "peak_demand": pval,
            "area": area,
        })
        row.append(area)
    winner_grid.append(row)

# Winner frequency table
print("\n=== Winner Frequency ===\n")
area_counts = Counter(r["area"] for r in results)
total = len(results)
for area, count in area_counts.most_common():
    print(f"  {area:<25} {count:>3}/{total}  ({100*count/total:.1f}%)")

# Stability stats
peak_lats = [r["peak_lat"] for r in results]
peak_lons = [r["peak_lon"] for r in results]
print(f"\nPeak latitude:  mean={np.mean(peak_lats):.4f}°N, "
      f"std={np.std(peak_lats):.4f}°")
print(f"Peak longitude: mean={np.mean(peak_lons):.4f}°E, "
      f"std={np.std(peak_lons):.4f}°")


# ============================================================
# 8. FIGURES
# ============================================================

# Plotting constants
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10

MODE_COLORS = {
    "tube": "#0019A8", "dlr": "#00AFAD",
    "elizabeth": "#6950A1", "rail": "#E87200"
}
MODE_MARKERS = {"tube": "s", "dlr": "D", "elizabeth": "^", "rail": "o"}
ASPECT = 111_000 / 69_000  # correct lat/lon aspect at ~51.5°N


# --- FIGURE 1: 2D Demand Heatmap ---
print("\nGenerating Figure 1: Demand heatmap...")

fig1, ax1 = plt.subplots(figsize=(10, 6))
im1 = ax1.pcolormesh(LON, LAT, demand_baseline, cmap='RdYlBu_r', shading='gouraud')
ax1.contour(LON, LAT, demand_baseline, levels=8, colors='k',
            linewidths=0.4, alpha=0.3)

# Plot stations by mode
for slat, slon, w, name, mode in STATIONS:
    if LAT_MIN <= slat <= LAT_MAX and LON_MIN <= slon <= LON_MAX:
        ax1.plot(slon, slat, MODE_MARKERS[mode], color=MODE_COLORS[mode],
                 markersize=7, markeredgecolor='white', markeredgewidth=0.8,
                 zorder=5)
        if w >= 3.5 or name in ["Kidbrooke", "Blackheath"]:
            ofs = (4, -10) if name in ["Woolwich Arsenal", "North Greenwich",
                                        "Charlton"] else (4, 4)
            ax1.annotate(name, (slon, slat), fontsize=6.5, color='#333',
                         xytext=ofs, textcoords='offset points',
                         fontfamily='serif', fontstyle='italic')

# Peak marker
ax1.plot(peak_lon, peak_lat, '*', color='#FF0000', markersize=14,
         markeredgecolor='white', markeredgewidth=1, zorder=10)
ax1.annotate(f'Peak demand\n({peak_lat:.3f}°N, {peak_lon:.3f}°E)',
             (peak_lon, peak_lat), fontsize=7.5, fontweight='bold',
             color='#CC0000', xytext=(12, -8), textcoords='offset points',
             arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1),
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='#CC0000', alpha=0.9))

# Proposed station
ax1.plot(0.005, 51.472, 'P', color='#00AA00', markersize=12,
         markeredgecolor='white', markeredgewidth=1, zorder=10)
ax1.annotate('Proposed station\n(Greenwich TC)',
             (0.005, 51.472), fontsize=7, fontweight='bold', color='#007700',
             xytext=(-65, 12), textcoords='offset points',
             arrowprops=dict(arrowstyle='->', color='#007700', lw=1),
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='#007700', alpha=0.9))

cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.85, pad=0.02)
cbar1.set_label('Unmet Transport Demand $D(x, y)$', fontsize=9)

legend_elements = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor='#0019A8',
           markersize=7, label='Underground'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='#6950A1',
           markersize=7, label='Elizabeth Line'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#00AFAD',
           markersize=7, label='DLR'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E87200',
           markersize=7, label='National Rail'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#FF0000',
           markersize=10, label='Peak demand'),
    Line2D([0], [0], marker='P', color='w', markerfacecolor='#00AA00',
           markersize=9, label='Proposed station'),
]
ax1.legend(handles=legend_elements, loc='lower left', fontsize=7,
           framealpha=0.9, edgecolor='#ccc', fancybox=False)

ax1.set_xlabel('Longitude (°E)')
ax1.set_ylabel('Latitude (°N)')
ax1.set_title('Unmet Transport Demand Surface — Royal Borough of Greenwich\n'
              r'$D(x,y) = \rho(x,y) \cdot (1 - A(x,y)/A_{\max})$, '
              r'baseline parameters ($\lambda = 800$m)',
              fontsize=11)
ax1.set_aspect(ASPECT)
plt.tight_layout()
plt.savefig('fig_demand_heatmap.png', dpi=300, bbox_inches='tight',
            facecolor='white')
plt.close()
print("  → Saved fig_demand_heatmap.png")


# --- FIGURE 2: 3D Surface ---
print("Generating Figure 2: 3D surface...")

fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection='3d')
surf = ax2.plot_surface(LON, LAT, demand_baseline, cmap='RdYlBu_r',
                        alpha=0.9, linewidth=0, antialiased=True,
                        rstride=1, cstride=1)

for slat, slon, w, name, mode in STATIONS:
    if LAT_MIN <= slat <= LAT_MAX and LON_MIN <= slon <= LON_MAX:
        li = np.argmin(np.abs(lats - slat))
        lj = np.argmin(np.abs(lons - slon))
        z = demand_baseline[li, lj] + 0.03
        ax2.scatter([slon], [slat], [z], c=MODE_COLORS[mode], s=25,
                    marker=MODE_MARKERS[mode], edgecolors='white',
                    linewidths=0.5, zorder=5)

ax2.scatter([peak_lon], [peak_lat], [demand_baseline.max() + 0.04],
            c='red', s=80, marker='*', edgecolors='white', linewidths=0.8,
            zorder=10)

ax2.set_xlabel('Longitude (°E)', fontsize=8, labelpad=8)
ax2.set_ylabel('Latitude (°N)', fontsize=8, labelpad=8)
ax2.set_zlabel('Demand $D(x,y)$', fontsize=8, labelpad=5)
ax2.set_title('3D Demand Surface — Greenwich Borough\n'
              r'Baseline parameters ($\lambda = 800$m, '
              'capacity-frequency weights)',
              fontsize=11, pad=15)
ax2.view_init(elev=32, azim=-55)
fig2.colorbar(surf, ax=ax2, shrink=0.55, pad=0.1, label='$D(x,y)$')
plt.tight_layout()
plt.savefig('fig_demand_3d.png', dpi=300, bbox_inches='tight',
            facecolor='white')
plt.close()
print("  → Saved fig_demand_3d.png")


# --- FIGURE 3: Sensitivity Grid ---
print("Generating Figure 3: Sensitivity grid...")

# Map area names to short codes for colouring
AREA_COLORS = {
    "Kidbrooke-Eltham": "#CC3333",
    "Thamesmead": "#3366CC",
    "Eltham": "#FF9900",
    "Kidbrooke corridor": "#CC3333",
    "Blackheath-Kidbrooke": "#CC3333",
    "Greenwich TC": "#33AA33",
}
SCHEME_NAMES = ["Equal", "Mild", "Baseline", "Steep"]

fig3, ax3 = plt.subplots(figsize=(8, 5))

for i, lam_val in enumerate(LAMBDAS):
    for j, scheme in enumerate(SCHEME_NAMES):
        area = winner_grid[i][j]
        color = AREA_COLORS.get(area, "#999999")
        rect = plt.Rectangle((j - 0.4, i - 0.35), 0.8, 0.7,
                              facecolor=color, edgecolor='white', linewidth=2)
        ax3.add_patch(rect)
        # Shorten long area names
        label = area.replace("Kidbrooke-Eltham", "Kidbrooke-\nEltham")
        ax3.text(j, i, label, ha='center', va='center',
                 fontsize=7, color='white', fontweight='bold')

        # Star the chosen baseline
        if lam_val == 800 and scheme == "Baseline":
            ax3.plot(j, i + 0.25, '*', color='yellow', markersize=12)

ax3.set_xticks(range(len(SCHEME_NAMES)))
ax3.set_xticklabels(SCHEME_NAMES, fontsize=9)
ax3.set_yticks(range(len(LAMBDAS)))
ax3.set_yticklabels([f'$\\lambda$ = {l}m' for l in LAMBDAS], fontsize=9)
ax3.set_xlim(-0.5, len(SCHEME_NAMES) - 0.5)
ax3.set_ylim(-0.5, len(LAMBDAS) - 0.5)
ax3.set_xlabel('Station Weight Scheme', fontsize=10)
ax3.set_ylabel('Catchment Decay Parameter', fontsize=10)
ax3.set_title('Sensitivity Analysis — Winning Area by Parameter Combination\n'
              '(★ = chosen baseline; '
              f'{area_counts.most_common(1)[0][1]}/20 favour '
              f'{area_counts.most_common(1)[0][0]})',
              fontsize=11)

legend_elements = [
    Patch(facecolor='#CC3333', label=f'Kidbrooke-Eltham ({area_counts.get("Kidbrooke-Eltham",0)+area_counts.get("Kidbrooke corridor",0)+area_counts.get("Blackheath-Kidbrooke",0)}/{total})'),
    Patch(facecolor='#3366CC', label=f'Thamesmead ({area_counts.get("Thamesmead",0)}/{total})'),
    Patch(facecolor='#FF9900', label=f'Eltham ({area_counts.get("Eltham",0)}/{total})'),
]
ax3.legend(handles=legend_elements, loc='upper right', fontsize=8,
           framealpha=0.95)
ax3.invert_yaxis()
plt.tight_layout()
plt.savefig('fig_sensitivity.png', dpi=300, bbox_inches='tight',
            facecolor='white')
plt.close()
print("  → Saved fig_sensitivity.png")


# --- FIGURE 4: Accessibility Field ---
print("Generating Figure 4: Accessibility field...")

_, acc_norm_plot = build_accessibility_field(lam=800)

fig4, ax4 = plt.subplots(figsize=(10, 6))
im4 = ax4.pcolormesh(LON, LAT, acc_norm_plot, cmap='viridis',
                      shading='gouraud')

for slat, slon, w, name, mode in STATIONS:
    if LAT_MIN <= slat <= LAT_MAX and LON_MIN <= slon <= LON_MAX:
        ax4.plot(slon, slat, MODE_MARKERS[mode], color='white',
                 markersize=7, markeredgecolor='#333',
                 markeredgewidth=0.8, zorder=5)
        if w >= 4.0:
            ax4.annotate(name, (slon, slat), fontsize=6.5, color='white',
                         xytext=(4, 4), textcoords='offset points',
                         fontstyle='italic',
                         path_effects=[pe.withStroke(linewidth=2,
                                                     foreground='black')])

ax4.annotate('Accessibility\ntrough', (0.038, 51.458), fontsize=9,
             color='#FF6666', fontweight='bold', ha='center',
             path_effects=[pe.withStroke(linewidth=2, foreground='black')])

cbar4 = plt.colorbar(im4, ax=ax4, shrink=0.85, pad=0.02)
cbar4.set_label(r'Normalised Accessibility $A(x,y) / A_{\max}$', fontsize=9)
ax4.set_xlabel('Longitude (°E)')
ax4.set_ylabel('Latitude (°N)')
ax4.set_title('Transport Accessibility Field — Royal Borough of Greenwich\n'
              r'Weighted sum of station influences '
              r'($\lambda = 800$m, capacity-frequency weights)',
              fontsize=11)
ax4.set_aspect(ASPECT)
plt.tight_layout()
plt.savefig('fig_accessibility.png', dpi=300, bbox_inches='tight',
            facecolor='white')
plt.close()
print("  → Saved fig_accessibility.png")


print("\n=== Done. 4 figures saved to working directory. ===")