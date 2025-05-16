import numpy as np
import matplotlib.pyplot as plt
import pykep as pk
from pykep.planet import jpl_lp
from pykep.core import epoch, DAY2SEC, MU_SUN, MU_EARTH
from pykep.core import lambert_problem
from pykep.trajopt._lambert import lambert_problem_multirev
from datetime import datetime, timedelta
from tqdm import tqdm

def mjd2000_to_datetime(mjd2000):
    base = datetime(2000, 1, 1, 12)
    return base + timedelta(days=float(mjd2000))

def datetime_to_mjd2000(dt):
    base = datetime(2000, 1, 1, 12)
    return (dt - base).total_seconds() / 86400.0

# Define planets
earth = jpl_lp("earth")
neptune = jpl_lp("neptune")

# Time grid
start_date = datetime(2033, 1, 1)
end_date = datetime(2038, 1, 1)
t0_start = datetime_to_mjd2000(start_date)
t0_end = datetime_to_mjd2000(end_date)

# Define grid steps
t0_grid = np.arange(t0_start, t0_end, 1)     # Departure dates
tof_grid = np.arange(1000, 4000, 1)          # Time of flight in days

# Prepare meshgrid for plotting
T0, TOF = np.meshgrid(t0_grid, tof_grid)
DV = np.full_like(T0, np.nan, dtype=float)

# Constants for Earth
r_leo = 6678e3       # m (Earth radius + 300 km)
vcirc = np.sqrt(MU_EARTH / r_leo)
vesc_sqr = 2 * MU_EARTH / r_leo

# Constants for Neptune
MU_NEPTUNE = 6.836529e15  # m^3/s^2
R_NEPTUNE = 24622e3        # m (mean radius)
r_neptune_orbit = R_NEPTUNE + 40000e3  # m
v_circ_neptune = np.sqrt(MU_NEPTUNE / r_neptune_orbit)
vesc_sqr_neptune = 2 * MU_NEPTUNE / r_neptune_orbit

# Main loop to compute porkchop data
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

            # Solve Lambert
            lambert = lambert_problem(r1, r2, tof * DAY2SEC, MU_SUN, cw=False, max_revs=0)
            lprob = lambert_problem_multirev(v1, lambert)

            v_dep = lprob.get_v1()[0]
            v_arr = lprob.get_v2()[0]

            vinf_dep = np.linalg.norm(np.array(v_dep) - np.array(v1))
            vinf_arr = np.linalg.norm(np.array(v_arr) - np.array(v2))

            # delta v at earth
            vlaunch = np.sqrt(vinf_dep**2 + vesc_sqr)
            delta_v_depart = vlaunch - vcirc

            # delta v at neptune
            v_arrival_hyp = np.sqrt(vinf_arr**2 + vesc_sqr_neptune)
            delta_v_arrival = v_arrival_hyp - v_circ_neptune

            # total_dv = vinf_dep + vinf_arr
            DV[i,j] = (delta_v_depart+delta_v_arrival) / 1000  # Convert to km/s

            # total_dv = vinf_dep + vinf_arr
            # DV[i, j] = total_dv / 1000  # km/s

        except:
            continue

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
levels = np.linspace(np.nanmin(DV), np.nanmax(DV), 30)
cs = ax.contourf(T0, TOF, DV, levels=levels, cmap='viridis')
cbar = plt.colorbar(cs)
cbar.set_label('Total Δv (km/s)')

# Format axes
t0_dates = [mjd2000_to_datetime(mjd).strftime('%Y-%m') for mjd in t0_grid]
ax.set_xticks(t0_grid[::4])
ax.set_xticklabels(t0_dates[::4], rotation=45)
ax.set_xlabel('Departure Date')
ax.set_ylabel('Time of Flight (days)')
ax.set_title('Earth → Neptune Porkchop Plot (Total Δv)')

plt.tight_layout()
plt.show()
