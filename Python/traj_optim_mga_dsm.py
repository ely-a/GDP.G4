import pykep as pk
import pygmo as pg
import numpy as np
from datetime import datetime
from pykep.planet import jpl_lp
from pykep.trajopt import mga_1dsm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting

def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days

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
    print(v_inf_entry)
    return v_inf_entry - v_para_neptune

def main():
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('neptune')]
    # seq = [jpl_lp('earth'), jpl_lp('neptune')]
    # seq = [jpl_lp('earth'), jpl_lp('mars'), jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('saturn'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2032, 3, 23)), to_mjd2000(datetime(2032, 3, 23))]
    tof = [
        [400, 466],  
        [2000, 2527],
    ]
    vinf = [10, 11]

    udp = mga_1dsm(
        seq=seq,
        t0=t0,
        tof=tof,
        vinf=vinf,
        add_vinf_dep=False,
        add_vinf_arr=True,
        tof_encoding="direct",     
        multi_objective=False
    )

    prob = pg.problem(udp)
    uda = pg.sade(gen=200)
    no_islands = 16
    archi = pg.archipelago(algo=uda, prob=prob, n=no_islands, pop_size=40)

    print(f"Running optimization on {no_islands} islands...")
    archi.evolve(50)
    archi.wait()

    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    x_best = archi.get_champions_x()[idx]

    # Compute v-infinity at departure and arrival
    DV_vec, _, _, _, _ = udp._compute_dvs(x_best)
    vinf_dep = np.linalg.norm([x_best[1], x_best[2], x_best[3]]) / 1000.0  # m/s to km/s
    vinf_arr = DV_vec[-1] / 1000.0  # m/s to km/s
    print(vinf_dep, vinf_arr)

    dv1 = delta_v_escape(vinf_dep)
    dv2 = delta_v_capture(vinf_arr)
    dsm_dv = float(sols[idx]) / 1000.0  # m/s to km/s

    total_dv = dsm_dv + dv1 + dv2 - vinf_arr

    udp.pretty(x_best)
    print(f"Best delta-v (DSM only): {dsm_dv:.3f} km/s")
    print(f"LEO escape delta-v: {dv1:.3f} km/s")
    print(f"Neptune capture delta-v: {dv2:.3f} km/s")
    print(f"Total mission delta-v: {total_dv:.3f} km/s")
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # udp.plot(x_best, ax=ax)
    # plt.show()

    _, _, _, ballistic_legs, ballistic_ep = udp._compute_dvs(x_best)

    eph_func = udp.get_eph_function(x_best)
    departure_epoch = x_best[0]
    tofs = udp._decode_tofs(x_best)
    arrival_epoch = sum(tofs) + x_best[0]
    times = np.linspace(departure_epoch, arrival_epoch, round(arrival_epoch-departure_epoch))  # 500 points

    from pykep.core import epoch

    planet_names = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    planets = [jpl_lp(name) for name in planet_names]

    # Prepare to store all positions
    all_positions = []

    for t in times:
        row = [t]
        for pl in planets:
            r, _ = pl.eph(epoch(t))
            # Convert from meters to kilometers
            row.extend([ri / 1e3 for ri in r])
        # Spacecraft position (already in meters, convert to km)
        sc_r, _ = eph_func(t)
        row.extend([ri / 1e3 for ri in sc_r])
        all_positions.append(row)

    # Write to file
    with open("spacecraft_trajectory.txt", "w") as f:
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

if __name__ == '__main__':
    main()
