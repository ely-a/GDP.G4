import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from tqdm import tqdm
import pykep as pk
from pykep.planet import jpl_lp
from pykep.core import epoch, DAY2SEC, MU_SUN, MU_EARTH
from pykep.core import lambert_problem, ic2par, propagate_lagrangian
from pykep.trajopt._lambert import lambert_problem_multirev
from matplotlib.dates import DateFormatter, MonthLocator

# --- Helper functions ---
def mjd2000_to_datetime(mjd2000):
    base = datetime(2000, 1, 1, 12)
    return base + timedelta(days=float(mjd2000))

def datetime_to_mjd2000(dt):
    base = datetime(2000, 1, 1, 12)
    return (dt - base).total_seconds() / 86400.0

# --- Define planets ---
earth = jpl_lp("earth")
neptune = jpl_lp("neptune")

# --- Time grid ---
start_date = datetime(2031, 5, 1)
end_date = datetime(2031, 8, 1)
t0_grid = np.arange(datetime_to_mjd2000(start_date), datetime_to_mjd2000(end_date), 1)  # every 5 days
tof_grid = np.arange(8.5 * 365.25, 9.5 * 365.25, 1)  # TOF from 1000 to 4050 days

# --- Meshgrids ---
T0, TOF = np.meshgrid(t0_grid, tof_grid)
DV = np.full_like(T0, np.nan, dtype=float)
DV_depart = np.full_like(T0, np.nan, dtype=float)
DV_arrival = np.full_like(T0, np.nan, dtype=float)
ARRIVAL = np.full_like(T0, np.nan, dtype=float)

# --- Earth departure orbit constants ---
r_leo = 6678e3  # Low Earth Orbit radius
vcirc = np.sqrt(MU_EARTH / r_leo)
vesc_sqr = 2 * MU_EARTH / r_leo

# --- Neptune arrival orbit constants ---
MU_NEPTUNE = 6.836529e15
R_NEPTUNE = 24622e3
r_neptune_orbit = R_NEPTUNE + 10000e3
v_circ_neptune = np.sqrt(MU_NEPTUNE / r_neptune_orbit)
v_para_neptune = np.sqrt(2*MU_NEPTUNE / r_neptune_orbit)
vesc_sqr_neptune = 2 * MU_NEPTUNE / r_neptune_orbit

# --- Main loop to compute porkchop plot data ---
for i in tqdm(range(T0.shape[0])):
    for j in range(T0.shape[1]):
        try:
            t0_val = T0[i, j]
            tof_val = TOF[i, j]
            t0 = epoch(float(t0_val))
            tof = float(tof_val)

            r1, v1 = earth.eph(t0)
            arrival_epoch = epoch(t0.mjd2000 + tof)
            r2, v2 = neptune.eph(arrival_epoch)

            lambert = lambert_problem(r1, r2, tof * DAY2SEC, MU_SUN, cw=False, max_revs=0)
            lprob = lambert_problem_multirev(v1, lambert)

            v_dep = lprob.get_v1()[0]
            v_arr = lprob.get_v2()[0]

            vinf_dep = np.linalg.norm(np.array(v_dep) - np.array(v1))
            vinf_arr = np.linalg.norm(np.array(v_arr) - np.array(v2))

            vlaunch = np.sqrt(vinf_dep**2 + vesc_sqr)
            delta_v_depart = vlaunch - vcirc

            v_arrival_hyp = np.sqrt(vinf_arr**2 + vesc_sqr_neptune)
            # delta_v_arrival = v_arrival_hyp - v_circ_neptune
            delta_v_arrival = v_arrival_hyp - v_para_neptune

            DV[i, j] = (delta_v_depart + delta_v_arrival) / 1000.0  # km/s
            DV_depart[i, j] = delta_v_depart / 1000.0  # km/s
            DV_arrival[i, j] = delta_v_arrival / 1000.0  # km/s
            ARRIVAL[i, j] = t0.mjd2000 + tof

        except:
            continue

# --- Convert MJD2000 to datetime for plotting ---
T0_dt = np.vectorize(mjd2000_to_datetime)(T0)
ARRIVAL_dt = np.vectorize(mjd2000_to_datetime)(ARRIVAL)

# --- User constraints ---
max_depart_dv = 10.0      # km/s, departure Δv
max_arrival_dv = 4.5     # km/s, arrival Δv
max_total_dv = max_depart_dv + max_arrival_dv      # km/s, total Δv
tof_min = 0.0            # years, min TOF
tof_max = 9.0 + 3.0/365.25            # years, max TOF

# --- Combined mask for all constraints ---
TOF_years = TOF / 365.25
valid_mask = (
    (DV <= max_total_dv) &
    (DV_depart <= max_depart_dv) &
    (DV_arrival <= max_arrival_dv) &
    (TOF_years >= tof_min) &
    (TOF_years <= tof_max)
)

# --- Individual masks for shading on each plot ---
mask_total = (DV <= max_total_dv) & (TOF_years >= tof_min) & (TOF_years <= tof_max)
mask_depart = (DV_depart <= max_depart_dv) & (TOF_years >= tof_min) & (TOF_years <= tof_max)
mask_arrival = (DV_arrival <= max_arrival_dv) & (TOF_years >= tof_min) & (TOF_years <= tof_max)

# --- Plotting: Total Δv ---
# fig, ax = plt.subplots(figsize=(12, 6))
# levels = np.linspace(np.nanmin(DV), np.nanmax(DV), 40)
# cs = ax.contourf(T0_dt, TOF_years, DV, levels=levels, cmap='viridis')
# cbar = plt.colorbar(cs)
# cbar.set_label('Total Δv (km/s)', fontsize=17)
# cbar.ax.tick_params(labelsize=17)
# # highlight = np.zeros_like(DV, dtype=float)
# # highlight[mask_total] = 1.0
# # ax.contourf(T0_dt, ARRIVAL_dt, highlight, levels=[0.5, 1.5], colors='none', hatches=['////'], alpha=0)
# cs_lines = ax.contour(T0_dt, TOF_years, DV, levels=levels, colors='k', linewidths=0.5)
# # ax.clabel(cs_lines, inline=True, fontsize=7, fmt="%.1f")
# tof_levels = np.arange(np.floor(np.nanmin(TOF_years)), np.ceil(np.nanmax(TOF_years))+0.2, 0.2)
# # cs_tof = ax.contour(
# #     T0_dt, ARRIVAL_dt, TOF_years,
# #     levels=tof_levels,
# #     colors='red',
# #     linestyles='dotted',
# #     linewidths=2.0
# # )
# # ax.clabel(cs_tof, fmt="%.1f yr", colors='red', fontsize=8)
# ax.xaxis.set_major_locator(MonthLocator(interval=4))
# ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
# # ax.yaxis.set_major_locator(MonthLocator(interval=6))
# # ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax.set_xlabel('Departure Date', fontsize=18)
# ax.set_ylabel('Time of Flight (years)', fontsize=18)
# ax.tick_params(axis='both', labelsize=17)
# fig.autofmt_xdate()
# plt.tight_layout()
# plt.show()

# # --- Plotting: Departure Δv only ---
# fig, ax = plt.subplots(figsize=(12, 6))
# levels = np.linspace(np.nanmin(DV_depart), np.nanmax(DV_depart), 30)
# cs = ax.contourf(T0_dt, ARRIVAL_dt, DV_depart, levels=levels, cmap='Blues')
# cbar = plt.colorbar(cs)
# cbar.set_label('Departure Δv (km/s)')
# highlight = np.zeros_like(DV_depart, dtype=float)
# highlight[mask_depart] = 1.0
# ax.contourf(T0_dt, ARRIVAL_dt, highlight, levels=[0.5, 1.5], colors='none', hatches=['////'], alpha=0)
# cs_lines = ax.contour(T0_dt, ARRIVAL_dt, DV_depart, levels=levels, colors='k', linewidths=0.5)
# ax.clabel(cs_lines, inline=True, fontsize=7, fmt="%.1f")
# cs_tof = ax.contour(
#     T0_dt, ARRIVAL_dt, TOF_years,
#     levels=tof_levels,
#     colors='red',
#     linestyles='dotted',
#     linewidths=2.0
# )
# ax.clabel(cs_tof, fmt="%.1f yr", colors='red', fontsize=8)
# ax.xaxis.set_major_locator(MonthLocator())
# ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax.yaxis.set_major_locator(MonthLocator(interval=3))
# ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax.set_xlabel('Departure Date')
# ax.set_ylabel('Arrival Date')
# ax.set_title('Earth → Neptune Porkchop Plot\nArrival Date vs Departure Date (Departure Δv)')
# fig.autofmt_xdate()
# plt.tight_layout()
# plt.show()

# # --- Plotting: Arrival Δv only ---
# fig, ax = plt.subplots(figsize=(12, 6))
# levels = np.linspace(np.nanmin(DV_arrival), np.nanmax(DV_arrival), 30)
# cs = ax.contourf(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels, cmap='Oranges')
# cbar = plt.colorbar(cs)
# cbar.set_label('Arrival Δv (km/s)')
# highlight = np.zeros_like(DV_arrival, dtype=float)
# highlight[mask_arrival] = 1.0
# ax.contourf(T0_dt, ARRIVAL_dt, highlight, levels=[0.5, 1.5], colors='none', hatches=['////'], alpha=0)
# cs_lines = ax.contour(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels, colors='k', linewidths=0.75)
# ax.clabel(cs_lines, inline=True, fontsize=8, fmt="%.1f")
# cs_tof = ax.contour(
#     T0_dt, ARRIVAL_dt, TOF_years,
#     levels=tof_levels,
#     colors='blue',
#     linestyles='dotted',
#     linewidths=2.0
# )
# ax.clabel(cs_tof, fmt="%.1f yr", colors='blue', fontsize=8)
# ax.xaxis.set_major_locator(MonthLocator())
# ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax.yaxis.set_major_locator(MonthLocator(interval=6))
# ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax.set_xlabel('Departure Date')
# ax.set_ylabel('Arrival Date')
# ax.set_title('Earth → Neptune Porkchop Plot\nArrival Date vs Departure Date (Arrival Δv)')
# fig.autofmt_xdate()
# plt.tight_layout()
# plt.show()

# import matplotlib.gridspec as gridspec

# # === Create Subplots ===
# fig = plt.figure(figsize=(15, 12))
# gs = gridspec.GridSpec(1, 2)
target_arrival_date = datetime(2040, 6, 2)

# # --- Departure Δv subplot ---
# ax1 = fig.add_subplot(gs[0])
# levels_depart = np.linspace(np.nanmin(DV_depart), np.nanmax(DV_depart), 20)
# cs1 = ax1.contourf(T0_dt, ARRIVAL_dt, DV_depart, levels=levels_depart, cmap='Blues')
# cbar1 = plt.colorbar(cs1, ax=ax1)
# cbar1.set_label('Departure Δv (km/s)', fontsize=18)
# cbar1.ax.tick_params(labelsize=15)
# cs_lines1 = ax1.contour(T0_dt, ARRIVAL_dt, DV_depart, levels=levels_depart, colors='k', linewidths=0.5)
# ax1.clabel(cs_lines1, inline=True, fontsize=12, fmt="%.1f")
# ax1.xaxis.set_major_locator(MonthLocator(interval=3))
# ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax1.yaxis.set_major_locator(MonthLocator(interval=3))
# ax1.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax1.set_xlabel('Departure Date', fontsize=18)
# ax1.tick_params(axis='both', labelsize=17)
# ax1.set_ylabel('Arrival Date', fontsize=18)
# ax1.axhline(target_arrival_date, color='red', linestyle='--', linewidth=1.5, label='Target: 2040-06-02')
# ax1.legend(loc='upper left', fontsize=15)

# # --- Arrival Δv subplot ---
# ax2 = fig.add_subplot(gs[1])
# levels_arrival = np.linspace(np.nanmin(DV_arrival), np.nanmax(DV_arrival), 20)
# cs2 = ax2.contourf(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels_arrival, cmap='Blues')
# cbar2 = plt.colorbar(cs2, ax=ax2)
# cbar2.set_label('Arrival Δv (km/s)', fontsize=18)
# cbar2.ax.tick_params(labelsize=15)
# cs_lines2 = ax2.contour(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels_arrival, colors='k', linewidths=0.5)
# ax2.clabel(cs_lines2, inline=True, fontsize=12, fmt="%.1f")
# ax2.xaxis.set_major_locator(MonthLocator(interval=3))
# ax2.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax2.yaxis.set_major_locator(MonthLocator(interval=3))
# ax2.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
# ax2.set_xlabel('Departure Date', fontsize = 18)
# ax2.set_ylabel('Arrival Date', fontsize=18)
# ax2.tick_params(axis='both', labelsize=17)
# ax2.axhline(target_arrival_date, color='red', linestyle='--', linewidth=1.5, label='Target: 2040-06-02')
# ax2.legend(loc='upper left', fontsize=15)


# plt.tight_layout()
# plt.subplots_adjust(bottom=0.1)  # Manually pad the bottom
# plt.show()

# --- Find and print transfer windows for all constraints ---
flat_dep_dates = T0_dt[valid_mask].flatten()
flat_arr_dates = ARRIVAL_dt[valid_mask].flatten()
flat_total_dv = DV[valid_mask].flatten()
flat_depart_dv = DV_depart[valid_mask].flatten()
flat_arrival_dv = DV_arrival[valid_mask].flatten()
flat_tof = TOF[valid_mask].flatten()

# Sort by departure date
sort_idx = np.argsort(flat_dep_dates)
flat_dep_dates = flat_dep_dates[sort_idx]
flat_arr_dates = flat_arr_dates[sort_idx]
flat_total_dv = flat_total_dv[sort_idx]
flat_depart_dv = flat_depart_dv[sort_idx]
flat_arrival_dv = flat_arrival_dv[sort_idx]
flat_tof = flat_tof[sort_idx]

threshold_days = 10  # adjust as needed
windows = []
if flat_dep_dates.size > 0:
    start_dep = flat_dep_dates[0]
    end_dep = flat_dep_dates[0]
    start_arr = flat_arr_dates[0]
    end_arr = flat_arr_dates[0]
    indices = [0]
    for i in range(1, len(flat_dep_dates)):
        gap = (flat_dep_dates[i] - flat_dep_dates[i-1]).days
        if gap > threshold_days:
            windows.append((start_dep, end_dep, start_arr, end_arr, indices))
            start_dep = flat_dep_dates[i]
            end_dep = flat_dep_dates[i]
            start_arr = flat_arr_dates[i]
            end_arr = flat_arr_dates[i]
            indices = [i]
        else:
            end_dep = flat_dep_dates[i]
            end_arr = flat_arr_dates[i]
            indices.append(i)
    windows.append((start_dep, end_dep, start_arr, end_arr, indices))

    # Print all windows and optimal departure in each
    for idx, (dep_start, dep_end, arr_start, arr_end, indices) in enumerate(windows, 1):
        print(f"Transfer window {idx}:")
        print(f"  Departure: {dep_start.strftime('%Y-%m-%d')} to {dep_end.strftime('%Y-%m-%d')}")
        print(f"  Arrival:   {arr_start.strftime('%Y-%m-%d')} to {arr_end.strftime('%Y-%m-%d')}")
        if indices:
            window_total_dv = flat_total_dv[indices]
            window_depart_dv = flat_depart_dv[indices]
            window_arrival_dv = flat_arrival_dv[indices]
            window_dep_dates = flat_dep_dates[indices]
            window_arr_dates = flat_arr_dates[indices]
            window_tof = flat_tof[indices]
            min_idx = np.argmin(window_total_dv)
            best_dep = window_dep_dates[min_idx]
            best_arr = window_arr_dates[min_idx]
            best_total_dv = window_total_dv[min_idx]
            best_depart_dv = window_depart_dv[min_idx]
            best_arrival_dv = window_arrival_dv[min_idx]
            best_tof = window_tof[min_idx]
            print(f"  Optimal departure: {best_dep.strftime('%Y-%m-%d')}")
            print(f"    Arrival:         {best_arr.strftime('%Y-%m-%d')}")
            print(f"    TOF:             {best_tof:.1f} days ({best_tof/365.25:.2f} years)")
            print(f"    Total Δv:        {best_total_dv:.4f} km/s")
            print(f"    Departure Δv:    {best_depart_dv:.4f} km/s")
            print(f"    Arrival Δv:      {best_arrival_dv:.4f} km/s")

            # Find and print the minimum TOF trajectory in this window
            min_tof_idx = np.argmin(window_tof)
            min_tof_dep = window_dep_dates[min_tof_idx]
            min_tof_arr = window_arr_dates[min_tof_idx]
            min_tof_total_dv = window_total_dv[min_tof_idx]
            min_tof_depart_dv = window_depart_dv[min_tof_idx]
            min_tof_arrival_dv = window_arrival_dv[min_tof_idx]
            min_tof_val = window_tof[min_tof_idx]
            print(f"  Minimum TOF departure: {min_tof_dep.strftime('%Y-%m-%d')}")
            print(f"    Arrival:             {min_tof_arr.strftime('%Y-%m-%d')}")
            print(f"    TOF:                 {min_tof_val:.1f} days ({min_tof_val/365.25:.2f} years)")
            print(f"    Total Δv:            {min_tof_total_dv:.2f} km/s")
            print(f"    Departure Δv:        {min_tof_depart_dv:.2f} km/s")
            print(f"    Arrival Δv:          {min_tof_arrival_dv:.2f} km/s")
        else:
            print("  No valid departures in this window.")
else:
    print("No valid transfer found for the given constraints.")

# Only do this for the first window with valid solutions
if windows:
    dep_start, dep_end, arr_start, arr_end, indices = windows[0]
    window_total_dv = flat_total_dv[indices]
    window_dep_dates = flat_dep_dates[indices]
    window_arr_dates = flat_arr_dates[indices]
    min_idx = np.argmin(window_total_dv)
    best_dep = window_dep_dates[min_idx]
    best_arr = window_arr_dates[min_idx]

    # Get initial and final state for the optimal transfer
    t0 = epoch(datetime_to_mjd2000(best_dep))
    tof = (best_arr - best_dep).days
    tof_sec = tof * 86400.0

    r1, v1 = earth.eph(t0)
    arrival_epoch = epoch(t0.mjd2000 + tof)
    r2, v2 = neptune.eph(arrival_epoch)

    lambert = lambert_problem(r1, r2, tof_sec, MU_SUN, cw=False, max_revs=0)
    v1_lam = lambert.get_v1()[0]
    v2_lam = lambert.get_v2()[0]
    vinf_dep = np.linalg.norm(np.array(v1_lam) - np.array(v1))
    vinf_arr = np.linalg.norm(np.array(v2_lam) - np.array(v2))

    # Compute orbital elements at departure
    elements = ic2par(r1, v1_lam, MU_SUN)
    a, e, i, raan, argp, ta = elements

    # Compute specific angular momentum
    h_vec = np.cross(r1, v1_lam)
    h_mag = np.linalg.norm(h_vec)

print("\nFirst optimal transfer details:")
print(f"  Departure date: {best_dep.strftime('%Y-%m-%d')} (MJD2000: {datetime_to_mjd2000(best_dep):.2f})")
print(f"  Arrival date:   {best_arr.strftime('%Y-%m-%d')} (MJD2000: {datetime_to_mjd2000(best_arr):.2f})")
print("\nOrbital elements of the transfer orbit at departure (first optimal solution):")
print(f"  Semi-major axis (a):      {a/1e3:.3f} km")
print(f"  Eccentricity (e):         {e:.6f}")
print(f"  Inclination (i):          {np.degrees(i):.3f} deg")
print(f"  RAAN:                     {np.degrees(raan):.3f} deg")
print(f"  Argument of periapsis:    {np.degrees(argp):.3f} deg")
print(f"  True anomaly:             {np.degrees(ta):.3f} deg")
print(f"  Specific angular momentum (h): {h_mag/1e6:.6f} km²/s")
print(f"  v-infinity at Earth:   {vinf_dep/1e3:.3f} km/s")
print(f"  v-infinity at Neptune: {vinf_arr/1e3:.3f} km/s")
v_entry = np.sqrt(vinf_arr**2 + 2 * MU_NEPTUNE / r_neptune_orbit)
print(f"  v at Neptune (entry): {v_entry/1e3:.3f} km/s")

# Number of points along the trajectory
N_points = 2000

# Generate epochs between departure and arrival (in MJD2000)
dep_mjd2000 = datetime_to_mjd2000(best_dep)
arr_mjd2000 = datetime_to_mjd2000(best_arr)
traj_epochs = np.linspace(dep_mjd2000, arr_mjd2000, N_points+1)[1:]

# Propagate spacecraft positions
sc_positions = []
for t_mjd in traj_epochs:
    dt_sec = (t_mjd - dep_mjd2000) * DAY2SEC
    r_vec, _ = propagate_lagrangian(r1, v1_lam, dt_sec, MU_SUN)
    sc_positions.append([ri / 1e3 for ri in r_vec])  # convert to km

sc_positions = np.array(sc_positions)  # shape (N_points, 3)

# --- Export trajectory and planet positions to text file ---

planet_names = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
planets = [jpl_lp(name) for name in planet_names]

# Add the departure epoch and positions as the first row
row0 = [dep_mjd2000]
for pl in planets:
    r, _ = pl.eph(epoch(dep_mjd2000))
    row0.extend([ri / 1e3 for ri in r])  # meters to km
row0.extend([ri / 1e3 for ri in r1])  # spacecraft at departure, in km
all_positions = [row0]

for t, sc_pos in zip(traj_epochs, sc_positions):
    row = [t]
    for pl in planets:
        r, _ = pl.eph(epoch(t))
        row.extend([ri / 1e3 for ri in r])  # meters to km
    row.extend(sc_pos)  # already in km
    all_positions.append(row)

with open("spacecraft_trajectory_direct.txt", "w") as f:
    # Header
    header = "time_mjd2000"
    for name in planet_names:
        header += f" {name}_x {name}_y {name}_z"
    header += " sc_x sc_y sc_z\n"
    f.write(header)
    # Data
    for row in all_positions:
        f.write(" ".join(f"{val:.8e}" for val in row) + "\n")

print("Exported full trajectory and planet positions for MATLAB.")