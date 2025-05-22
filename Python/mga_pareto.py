import pykep as pk
import pygmo as pg
from pykep.planet import jpl_lp
from pykep.trajopt import mga
from datetime import datetime, timedelta
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days

# Constants for Earth and Neptune
r_earth = 6378.0          # km
mu_earth = 398600.4418    # km^3/s^2
leo_altitude = 300.0      # km
r_leo = r_earth + leo_altitude

mu_neptune = 6835100      # km^3/s^2
r_neptune = 24622.0       # km
neptune_orbit_alt = 10000.0  # km
r_neptune_orbit = r_neptune + neptune_orbit_alt

# Δv to escape LEO (300 km) into hyperbolic trajectory
def delta_v_escape(vinf):  # vinf in km/s
    v_circ = np.sqrt(mu_earth / r_leo)
    v_peri_hyp = np.sqrt(vinf**2 + 2 * mu_earth / r_leo)
    return v_peri_hyp - v_circ

# Δv to capture into circular orbit around Neptune
def delta_v_capture(vinf):  # vinf in km/s
    # v_circ = np.sqrt(mu_neptune / r_neptune_orbit)
    v_para_neptune = np.sqrt(2*mu_neptune / r_neptune_orbit)
    v_inf_entry = np.sqrt(vinf**2 + 2 * mu_neptune / r_neptune_orbit)
    return v_inf_entry - v_para_neptune

def mjd2000_to_datetime(mjd2000):
    return datetime(2000, 1, 1, 12) + timedelta(days=float(mjd2000))

# ------------------------------------
# Define MGA Optimization Problem
# ------------------------------------
seq = [jpl_lp('earth'), jpl_lp('mars'), jpl_lp('jupiter'), jpl_lp('neptune')]
N_legs = len(seq) - 1
t0 = [to_mjd2000(datetime(2029, 5, 13)), to_mjd2000(datetime(2032, 5, 13))]
tof = [[100, 200], [300, 600], [1000, 3500]]
vinf = 7.0

udp = mga(
    seq=seq,
    t0=t0,
    tof=tof,
    vinf=vinf,
    tof_encoding="direct",
    multi_objective=True
)

prob = pg.problem(udp)
algo = pg.algorithm(pg.nsga2(gen=50))

pop = pg.population(prob, size=200)
for _ in range(10):  # 500 generations total
    pop = algo.evolve(pop)

f_vals = pop.get_f()         # [dv_total (sum of impulses), TOF]
x_vals = pop.get_x()         # decision vectors
adjusted_dvs = []
dep_mjd2000 = []
arr_mjd2000 = []

for xi, fi in zip(x_vals, f_vals):
    # For mga, direct encoding:
    # [t0, T1, T2, ..., v_inf_dep_x, v_inf_dep_y, v_inf_dep_z, v_inf_arr_x, v_inf_arr_y, v_inf_arr_z]
    vinf_dep = np.linalg.norm(xi[1+N_legs:1+N_legs+3])  # m/s
    vinf_arr = np.linalg.norm(xi[-3:])                  # m/s

    dv1 = delta_v_escape(vinf_dep/1000.0)  # Convert m/s → km/s
    dv2 = delta_v_capture(vinf_arr/1000.0)  # Convert m/s → km/s
    total_dv = fi[0] / 1000.0 + dv1 + dv2  # fi[0] is total Δv in m/s

    adjusted_dvs.append(total_dv)
    dep_mjd2000.append(xi[0])
    tof_sum = sum(xi[1:1+N_legs])
    arr_mjd2000.append(xi[0] + tof_sum)

adjusted_dvs = np.array(adjusted_dvs)
# dep_mjd2000 = np.array(dep_mjd2000)
# arr_mjd2000 = np.array(arr_mjd2000)

# # 3D scatter: x=Δv, y=TOF, z=departure date (MJD2000)
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')
# sc = ax.scatter(adjusted_dvs, f_vals[:, 1], dep_mjd2000, c=adjusted_dvs, cmap='viridis', s=20)
# ax.set_xlabel("Total Δv incl. LEO injection and Neptune capture [km/s]")
# ax.set_ylabel("Time of Flight [days]")
# ax.set_zlabel("Departure Date (MJD2000)")
# z_ticks = np.linspace(dep_mjd2000.min(), dep_mjd2000.max(), 8)
# ax.set_zticks(z_ticks)
# ax.set_zticklabels([mjd2000_to_datetime(tick).strftime('%Y-%m') for tick in z_ticks])
# fig.colorbar(sc, ax=ax, label='Total Δv [km/s]')
# ax.set_title("Pareto Front: Δv vs TOF vs Departure Date (MGA)")
# plt.tight_layout()
# plt.show()

# # Δv vs hours since first departure
# first_dep = dep_mjd2000.min()
# hours_since_first = (dep_mjd2000 - first_dep) * 24  # 1 day = 24 hours
# plt.figure(figsize=(10, 6))
# plt.scatter(hours_since_first, adjusted_dvs, c=adjusted_dvs, cmap='viridis', s=20)
# plt.xlabel("Hours since first departure date")
# plt.ylabel("Total Δv incl. LEO injection and Neptune capture [km/s]")
# plt.title("Δv vs Hours Since First Departure (MGA)")
# plt.colorbar(label='Total Δv [km/s]')
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# # Arrival vs Departure date, Δv as contour
# def get_tick_positions(mjd_min, mjd_max, months=6):
#     dt_min = mjd2000_to_datetime(mjd_min)
#     dt_max = mjd2000_to_datetime(mjd_max)
#     ticks = []
#     dt = datetime(dt_min.year, dt_min.month, 1, 12)
#     while dt <= dt_max:
#         ticks.append((dt - datetime(2000, 1, 1, 12)).days)
#         # increment by N months
#         m = dt.month - 1 + months
#         y = dt.year + m // 12
#         m = m % 12 + 1
#         dt = datetime(y, m, 1, 12)
#     return ticks

# fig, ax = plt.subplots(figsize=(10, 6))
# levels = np.linspace(np.nanmin(adjusted_dvs), np.nanmax(adjusted_dvs), 30)
# cs = ax.tricontourf(dep_mjd2000, arr_mjd2000, adjusted_dvs, levels=levels, cmap='viridis')
# cbar = plt.colorbar(cs)
# cbar.set_label("Total Δv incl. LEO injection and Neptune capture [km/s]")

# x_ticks = get_tick_positions(dep_mjd2000.min(), dep_mjd2000.max(), months=6)
# y_ticks = get_tick_positions(arr_mjd2000.min(), arr_mjd2000.max(), months=6)
# ax.set_xticks(x_ticks)
# ax.set_yticks(y_ticks)
# ax.set_xticklabels([mjd2000_to_datetime(tick).strftime('%Y-%m') for tick in x_ticks], rotation=45)
# ax.set_yticklabels([mjd2000_to_datetime(tick).strftime('%Y-%m') for tick in y_ticks])
# ax.set_xlabel("Departure Date")
# ax.set_ylabel("Arrival Date")
# ax.set_title("Pareto MGA: Arrival vs Departure Date (Δv as contour)")
# plt.tight_layout()
# plt.show()

best_idx = np.argmin(adjusted_dvs)
best_x = x_vals[best_idx]

print("Best Solution (detailed):")
udp.pretty(best_x)