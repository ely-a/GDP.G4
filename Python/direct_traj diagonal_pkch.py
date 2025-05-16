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
t0_grid = np.arange(datetime_to_mjd2000(start_date), datetime_to_mjd2000(end_date), 1)  # every 5 days
tof_grid = np.arange(2500, 4050, 1)  # TOF from 1000 to 4050 days

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

# --- Plotting: Total Δv ---
fig, ax = plt.subplots(figsize=(12, 6))
levels = np.linspace(np.nanmin(DV), np.nanmax(DV), 40)
cs = ax.contourf(T0_dt, ARRIVAL_dt, DV, levels=levels, cmap='viridis')
cbar = plt.colorbar(cs)
cbar.set_label('Total Δv (km/s)')
cs_lines = ax.contour(T0_dt, ARRIVAL_dt, DV, levels=levels, colors='k', linewidths=0.5)
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

ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.yaxis.set_major_locator(MonthLocator(interval=6))
ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.set_xlabel('Departure Date')
ax.set_ylabel('Arrival Date')
ax.set_title('Earth → Neptune Porkchop Plot\nArrival Date vs Departure Date (Total Δv)')
fig.autofmt_xdate()
plt.tight_layout()
plt.show()

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

# --- Plotting: Arrival Δv only ---
fig, ax = plt.subplots(figsize=(12, 6))
levels = np.linspace(np.nanmin(DV_arrival), np.nanmax(DV_arrival), 30)
cs = ax.contourf(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels, cmap='Oranges')
cbar = plt.colorbar(cs)
cbar.set_label('Arrival Δv (km/s)')
cs_lines = ax.contour(T0_dt, ARRIVAL_dt, DV_arrival, levels=levels, colors='k', linewidths=0.75)

# Add TOF (years) red dotted contours
TOF_years = TOF / 365.25
tof_levels = np.arange(np.floor(np.nanmin(TOF_years)), np.ceil(np.nanmax(TOF_years))+0.5, 0.5)
cs_tof = ax.contour(
    T0_dt, ARRIVAL_dt, TOF_years,
    levels=tof_levels,
    colors='blue',
    linestyles='dotted',
    linewidths=2.0
)
ax.clabel(cs_tof, fmt="%.1f yr", colors='blue', fontsize=8)

ax.clabel(cs_lines, inline=True, fontsize=8, fmt="%.1f")
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.yaxis.set_major_locator(MonthLocator(interval=6))
ax.yaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.set_xlabel('Departure Date')
ax.set_ylabel('Arrival Date')
ax.set_title('Earth → Neptune Porkchop Plot\nArrival Date vs Departure Date (Arrival Δv)')
fig.autofmt_xdate()
plt.tight_layout()
plt.show()
