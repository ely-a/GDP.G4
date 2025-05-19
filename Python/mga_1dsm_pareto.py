import pykep as pk
import pygmo as pg
from pykep.planet import jpl_lp
from pykep.trajopt import mga_1dsm
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
neptune_orbit_alt = 40000.0  # km
r_neptune_orbit = r_neptune + neptune_orbit_alt

# Δv to escape LEO (300 km) into hyperbolic trajectory
def delta_v_escape(vinf):  # vinf in km/s
    v_circ = np.sqrt(mu_earth / r_leo)
    v_peri_hyp = np.sqrt(vinf**2 + 2 * mu_earth / r_leo)
    return v_peri_hyp - v_circ

# Δv to capture into circular orbit around Neptune
def delta_v_capture(vinf):  # vinf in km/s
    v_circ = np.sqrt(mu_neptune / r_neptune_orbit)
    v_inf_entry = np.sqrt(vinf**2 + 2 * mu_neptune / r_neptune_orbit)
    return v_inf_entry - v_circ

# ------------------------------------
# Define MGA-1DSM Optimization Problem
# ------------------------------------
seq = [jpl_lp('earth'), jpl_lp('mars'), jpl_lp('jupiter'), jpl_lp('neptune')]
t0 = [to_mjd2000(datetime(2029, 5, 13)), to_mjd2000(datetime(2032, 5, 13))]
tof = [[100, 200], [300, 600], [1000, 3500]]
vinf = [6, 9.5]

udp = mga_1dsm(
    seq=seq,
    t0=t0,
    tof=tof,
    vinf=vinf,
    add_vinf_dep=False,   # include v∞ in Δv at departure
    add_vinf_arr=False,   # include v∞ in Δv at arrival
    tof_encoding="direct",
    multi_objective=True
)

prob = pg.problem(udp)
algo = pg.algorithm(pg.nsga2(gen=50))

pop = pg.population(prob, size=200)
for _ in range(10):  # 500 generations total
    pop = algo.evolve(pop)

f_vals = pop.get_f()         # [dv_total (DSM+v∞), TOF]
x_vals = pop.get_x()         # decision vectors
adjusted_dvs = []

# --------------------------------------------
# Adjust Δv: Add LEO → v∞ + v∞ → capture burns
# --------------------------------------------
for xi, fi in zip(x_vals, f_vals):
    DV, _, _, _, _ = udp._compute_dvs(xi)
    vinf_dep = np.linalg.norm([xi[1], xi[2], xi[3]])  # m/s
    vinf_arr = DV[-1]  # m/s

    dv1 = delta_v_escape(vinf_dep/1000.0)  # Convert m/s → km/s
    dv2 = delta_v_capture(vinf_arr/1000.0)  # Convert m/s → km/s
    dsm_dv = fi[0] / 1000.0              # Convert m/s → km/s

    total_dv = dsm_dv + dv1 + dv2
    adjusted_dvs.append(total_dv)

adjusted_dvs = np.array(adjusted_dvs)

# dep_mjd2000 = np.array([xi[0] for xi in x_vals])

# def mjd2000_to_datetime(mjd2000):
#     return datetime(2000, 1, 1, 12) + timedelta(days=float(mjd2000))
# dep_dates = np.array([mjd2000_to_datetime(mjd) for mjd in dep_mjd2000])

# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')

# # 3D scatter: x=Δv, y=TOF, z=departure date (MJD2000)
# sc = ax.scatter(adjusted_dvs, f_vals[:, 1], dep_mjd2000, c=adjusted_dvs, cmap='viridis', s=20)

# ax.set_xlabel("Total Δv incl. LEO injection and Neptune capture [km/s]")
# ax.set_ylabel("Time of Flight [days]")
# ax.set_zlabel("Departure Date (MJD2000)")

# # Optional: set z-ticks as dates
# z_ticks = np.linspace(dep_mjd2000.min(), dep_mjd2000.max(), 8)
# ax.set_zticks(z_ticks)
# ax.set_zticklabels([mjd2000_to_datetime(tick).strftime('%Y-%m') for tick in z_ticks])

# fig.colorbar(sc, ax=ax, label='Total Δv [km/s]')
# ax.set_title("Pareto Front: Δv vs TOF vs Departure Date (MGA-1DSM)")
# plt.tight_layout()
# plt.show()
# ---------------------
# Plot adjusted results
# ---------------------
plt.figure(figsize=(10, 6))
plt.scatter(adjusted_dvs, f_vals[:, 1], s=15, c='green')
plt.xlabel("Total Δv incl. LEO injection and Neptune capture [km/s]")
plt.ylabel("Time of Flight [days]")
plt.title("Adjusted Pareto Front: Earth → Mars → Jupiter → Neptune (MGA-1DSM)")
plt.grid(True)
plt.tight_layout()
plt.show()

# first_dep = dep_mjd2000.min()
# hours_since_first = (dep_mjd2000 - first_dep) * 24  # 1 day = 24 hours

# plt.figure(figsize=(10, 6))
# plt.scatter(hours_since_first, adjusted_dvs, c=adjusted_dvs, cmap='viridis', s=20)
# plt.xlabel("Hours since first departure date")
# plt.ylabel("Total Δv incl. LEO injection and Neptune capture [km/s]")
# plt.title("Δv vs Hours Since First Departure (MGA-1DSM)")
# plt.colorbar(label='Total Δv [km/s]')
# plt.grid(True)
# plt.tight_layout()
# plt.show()