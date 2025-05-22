import pykep as pk
import pygmo as pg
import numpy as np
from datetime import datetime
from pykep.planet import jpl_lp
from pykep.trajopt import mga

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
    return v_inf_entry - v_para_neptune

def main():
    # Define planets and bounds
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2032, 3, 23)), to_mjd2000(datetime(2032, 3, 24))]  # launch window
    tof = [
        [400, 500],   # Earth->Jupiter leg (days)
        [2000, 2527], # Jupiter->Neptune leg (days)
    ]
    vinf = 11.0  # km/s (launch v-infinity)

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
    uda = pg.sade(gen=200)
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

if __name__ == "__main__":
    main()