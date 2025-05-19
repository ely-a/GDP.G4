import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from tqdm import tqdm
import pykep as pk
from pykep.planet import jpl_lp
from pykep.core import epoch, DAY2SEC, MU_SUN, MU_EARTH
from pykep.core import lambert_problem
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
start_date = datetime(2033, 1, 1)
end_date = datetime(2035, 1, 1)
t0_grid = np.arange(datetime_to_mjd2000(start_date), datetime_to_mjd2000(end_date), 1)  # every day
tof_grid = np.arange(2500, 4050, 1)  # TOF from 2500 to 4050 days

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
r_neptune_orbit = R_NEPTUNE + 40000e3
v_circ_neptune = np.sqrt(MU_NEPTUNE / r_neptune_orbit)
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
            delta_v_arrival = v_arrival_hyp - v_circ_neptune

            DV[i, j] = (delta_v_depart + delta_v_arrival) / 1000.0  # km/s
            DV_depart[i, j] = delta_v_depart / 1000.0  # km/s
            DV_arrival[i, j] = delta_v_arrival / 1000.0  # km/s
            ARRIVAL[i, j] = t0.mjd2000 + tof

        except:
            continue

# --- Convert MJD2000 to datetime for plotting ---
T0_dt = np.vectorize(mjd2000_to_datetime)(T0)
ARRIVAL_dt = np.vectorize(mjd2000_to_datetime)(ARRIVAL)

# --- Plotting: Departure Δv only ---
fig, ax = plt.subplots(figsize=(12, 6))
levels = np.linspace(np.nanmin(DV_depart), np.nanmax(DV_depart), 30)
cs = ax.contourf(T0_dt, ARRIVAL_dt, DV_depart, levels=levels, cmap='Blues')
cbar = plt.colorbar(cs)
cbar.set_label('Departure Δv (km/s)')
cs_lines = ax.contour(T0_dt, ARRIVAL_dt, DV_depart, levels=levels, colors='k', linewidths=0.5)
ax.clabel(cs_lines, inline=True, fontsize=7, fmt="%.1f")

# Add TOF (years) red dotted contours
TOF_years = TOF / 365.25
tof_levels = np.arange(np.floor(np.nanmin(TOF_years)), np.ceil(np.nanmax(TOF_years))+0.5, 0.5)
cs_tof = ax.contour(
    T0_dt, ARRIVAL_dt, TOF_years,
    levels=tof_levels,
    colors='red',
    linestyles='dotted',
    linewidths=2.0
)
ax.clabel(cs_tof, fmt="%.1f yr", colors='red', fontsize=8)

# --- Constrain by departure Δv ---
max_depart_dv = 9.5  # Example: set your max allowed departure delta-v in km/s
tof_min = 0.0        # years
tof_max = 11.0       # years

# Mask for solutions that meet both constraints
valid_mask = (DV_depart <= max_depart_dv) & (TOF_years >= tof_min) & (TOF_years <= tof_max)

# Highlight valid region on departure Δv plot
highlight = np.zeros_like(DV_depart, dtype=float)
highlight[valid_mask] = 1.0
ax.contourf(T0_dt, ARRIVAL_dt, highlight, levels=[0.5, 1.5], colors='none', hatches=['////'], alpha=0)

ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.yaxis.set_major_locator(MonthLocator(interval=6))
ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.set_xlabel('Departure Date')
ax.set_ylabel('Arrival Date')
ax.set_title('Earth → Neptune Porkchop Plot\nArrival Date vs Departure Date (Departure Δv)')
fig.autofmt_xdate()
plt.tight_layout()
plt.show()

# --- Find and print transfer windows for departure Δv constraint ---
flat_dep_dates = T0_dt[valid_mask].flatten()
flat_arr_dates = ARRIVAL_dt[valid_mask].flatten()

# Sort by departure date
sort_idx = np.argsort(flat_dep_dates)
flat_dep_dates = flat_dep_dates[sort_idx]
flat_arr_dates = flat_arr_dates[sort_idx]

# Group into windows (gap > threshold means new window)
threshold_days = 10  # adjust as needed
windows = []
if flat_dep_dates.size > 0:
    start_dep = flat_dep_dates[0]
    end_dep = flat_dep_dates[0]
    start_arr = flat_arr_dates[0]
    end_arr = flat_arr_dates[0]
    for i in range(1, len(flat_dep_dates)):
        gap = (flat_dep_dates[i] - flat_dep_dates[i-1]).days
        if gap > threshold_days:
            windows.append((start_dep, end_dep, start_arr, end_arr))
            start_dep = flat_dep_dates[i]
            end_dep = flat_dep_dates[i]
            start_arr = flat_arr_dates[i]
            end_arr = flat_arr_dates[i]
        else:
            end_dep = flat_dep_dates[i]
            end_arr = flat_arr_dates[i]
    # Add last window
    windows.append((start_dep, end_dep, start_arr, end_arr))

    # Print all windows and optimal departure in each
    for idx, (dep_start, dep_end, arr_start, arr_end) in enumerate(windows, 1):
        print(f"Transfer window {idx}:")
        print(f"  Departure: {dep_start.strftime('%Y-%m-%d')} to {dep_end.strftime('%Y-%m-%d')}")
        print(f"  Arrival:   {arr_start.strftime('%Y-%m-%d')} to {arr_end.strftime('%Y-%m-%d')}")
        # Find all departures in this window
        in_window = (flat_dep_dates >= dep_start) & (flat_dep_dates <= dep_end)
        if np.any(in_window):
            dvs_in_window = DV_depart[valid_mask].flatten()[sort_idx][in_window]
            dep_dates_in_window = flat_dep_dates[in_window]
            arr_dates_in_window = flat_arr_dates[in_window]
            min_idx = np.argmin(dvs_in_window)
            best_dep = dep_dates_in_window[min_idx]
            best_arr = arr_dates_in_window[min_idx]
            best_dv = dvs_in_window[min_idx]
            print(f"  Optimal departure: {best_dep.strftime('%Y-%m-%d')} (Δv = {best_dv:.2f} km/s)")
            print(f"  Optimal arrival:   {best_arr.strftime('%Y-%m-%d')}")
        else:
            print("  No valid departures in this window.")
else:
    print("No valid transfer found for the given constraints.")

