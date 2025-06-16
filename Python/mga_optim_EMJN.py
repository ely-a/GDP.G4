import pykep as pk
import pygmo as pg
import numpy as np
from datetime import datetime
from pykep.planet import jpl_lp
from pykep.trajopt import mga
from pykep.core import propagate_lagrangian, ic2par, DAY2SEC, epoch

def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days

def datetime_to_mjd2000(dt):
    base = datetime(2000, 1, 1, 12)
    return (dt - base).total_seconds() / 86400.0

def delta_v_escape(vinf):  # vinf in km/s
    r_earth = 6378.0
    mu_earth = 398600.4418
    leo_altitude = 300.0
    r_leo = r_earth + leo_altitude
    v_circ = np.sqrt(mu_earth / r_leo)
    v_peri_hyp = np.sqrt(vinf**2 + 2 * mu_earth / r_leo)
    return v_peri_hyp - v_circ

def delta_v_capture(vinf):  # vinf in km/s
    mu_neptune = 6835100
    r_neptune = 24622.0
    neptune_orbit_alt = 10000.0
    r_neptune_orbit = r_neptune + neptune_orbit_alt
    v_para_neptune = np.sqrt(2*mu_neptune / r_neptune_orbit)
    v_inf_entry = np.sqrt(vinf**2 + 2 * mu_neptune / r_neptune_orbit)
    return v_inf_entry - v_para_neptune

def main():
    # Define planets and bounds
    seq = [jpl_lp('earth'), jpl_lp('mars'), jpl_lp('jupiter'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2031, 2, 8)), to_mjd2000(datetime(2031, 2, 9))]  # launch window
    tof = [
        [90, 100],   # Earth->Jupiter leg (days)
        [500, 600], # Jupiter->Neptune leg (days)
        [3300, 3333]
    ]
    vinf = 8.68  # km/s (launch v-infinity)

    # Create the MGA UDP
    udp = mga(
        seq=seq,
        t0=t0,
        tof=tof,
        vinf=vinf,
        multi_objective=False,
        tof_encoding='direct'
    )

    prob = pg.problem(udp)
    uda = pg.sade(gen=150)
    archi = pg.archipelago(algo=uda, prob=prob, n=8, pop_size=30)

    print("Running optimization...")
    archi.evolve(20)
    archi.wait()

    sols = archi.get_champions_f()
    idx = np.argmin(sols)
    x_best = archi.get_champions_x()[idx]

    # Print trajectory summary
    udp.pretty(x_best)

    DVlaunch, DVfb, DVarrival, l, DVlaunch_tot, T, ballistic_legs, ballistic_ep = udp._compute_dvs(x_best)

    vinf_dep = DVlaunch_tot / 1000.0  # km/s
    vinf_arr = DVarrival / 1000.0     # km/s

    dv1 = delta_v_escape(vinf_dep)
    dv2 = delta_v_capture(vinf_arr)
    total_dv = dv1 + dv2

    print(f"LEO escape delta-v: {dv1:.3f} km/s")
    print(f"Neptune capture delta-v: {dv2:.3f} km/s")
    print(f"Total mission delta-v (LEO escape + Neptune capture): {total_dv:.3f} km/s")
    # --- Propagate each leg of the MGA trajectory ---

    N_points = 2000
    traj_epochs = []
    sc_positions = []

    # Leg 1: Earth to Jupiter
    dep_mjd2000 = x_best[0]
    tof1 = x_best[1]
    arr1_mjd2000 = dep_mjd2000 + tof1

    r1, v1 = seq[0].eph(epoch(dep_mjd2000))
    r2, v2 = seq[1].eph(epoch(arr1_mjd2000))
    dt1 = tof1 * DAY2SEC
    lambert1 = pk.lambert_problem(r1, r2, dt1, pk.MU_SUN)
    v1_lam = lambert1.get_v1()[0]

    N1 = int(N_points * (tof1 / sum(x_best[1:])))
    epochs1 = np.linspace(dep_mjd2000, arr1_mjd2000, N1, endpoint=False)
    for t_mjd in epochs1:
        dt_sec = (t_mjd - dep_mjd2000) * DAY2SEC
        r_vec, _ = propagate_lagrangian(r1, v1_lam, dt_sec, pk.MU_SUN)
        sc_positions.append([ri / 1e3 for ri in r_vec])
        traj_epochs.append(t_mjd)

    # Leg 2: Jupiter to Neptune
    # Leg 2: Mars to Jupiter
    tof2 = x_best[2]
    arr2_mjd2000 = arr1_mjd2000 + tof2

    r3, v3 = seq[2].eph(epoch(arr2_mjd2000))
    dt2 = tof2 * DAY2SEC
    lambert2 = pk.lambert_problem(r2, r3, dt2, pk.MU_SUN)
    v2_lam = lambert2.get_v1()[0]

    N2 = int(N_points * (tof2 / sum(x_best[1:])))
    epochs2 = np.linspace(arr1_mjd2000, arr2_mjd2000, N2, endpoint=False)
    for idx, t_mjd in enumerate(epochs2):
        dt_sec = (t_mjd - arr1_mjd2000) * DAY2SEC
        r_vec, _ = propagate_lagrangian(r2, v2_lam, dt_sec, pk.MU_SUN)
        sc_positions.append([ri / 1e3 for ri in r_vec])
        traj_epochs.append(t_mjd)

    # Leg 3: Jupiter to Neptune
    tof3 = x_best[3]
    arr3_mjd2000 = arr2_mjd2000 + tof3

    r4, v4 = seq[3].eph(epoch(arr3_mjd2000))
    dt3 = tof3 * DAY2SEC
    lambert3 = pk.lambert_problem(r3, r4, dt3, pk.MU_SUN)
    v3_lam = lambert3.get_v1()[0]

    N3 = N_points - N1 - N2 + 1  # balance the total number of points
    epochs3 = np.linspace(arr2_mjd2000, arr3_mjd2000, N3)
    for idx, t_mjd in enumerate(epochs3):
        if idx == 0:
            continue
        dt_sec = (t_mjd - arr2_mjd2000) * DAY2SEC
        r_vec, _ = propagate_lagrangian(r3, v3_lam, dt_sec, pk.MU_SUN)
        sc_positions.append([ri / 1e3 for ri in r_vec])
        traj_epochs.append(t_mjd)

    sc_positions = np.array(sc_positions)
    traj_epochs = np.array(traj_epochs)


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

    with open("spacecraft_trajectory_mga.txt", "w") as f:
        # Header
        header = "time_mjd2000"
        for name in planet_names:
            header += f" {name}_x {name}_y {name}_z"
        header += " sc_x sc_y sc_z\n"
        f.write(header)
        # Data
        for row in all_positions:
            f.write(" ".join(f"{val:.8e}" for val in row) + "\n")

    print("Exported MGA trajectory and planet positions for MATLAB.")

if __name__ == "__main__":
    main()