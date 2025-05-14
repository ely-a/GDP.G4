import pykep as pk
import pygmo as pg
from datetime import datetime
from pykep.planet import jpl_lp
from pykep.trajopt import mga_1dsm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting

def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days

def main():
    seq = [jpl_lp('earth'), jpl_lp('mars'), jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('saturn'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2028, 5, 13)), to_mjd2000(datetime(2031, 5, 13))]
    tof = [
        [100, 250],     # Earth → Mars
        [150, 500],     # Mars → Earth
        [300, 800],     # Earth → Jupiter
        [400, 1000],    # Jupiter → Saturn
        [600, 1750]     # Saturn → Neptune
    ]
    vinf = [3, 10]

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
    no_islands = 8
    archi = pg.archipelago(algo=uda, prob=prob, n=no_islands, pop_size=20)

    print(f"Running optimization on {no_islands} islands...")
    archi.evolve(20)
    archi.wait()

    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    x_best = archi.get_champions_x()[idx]

    print("Best delta-v:", sols[idx])
    udp.pretty(x_best)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # udp.plot(x_best, ax=ax)
    # plt.show()

    print(udp._decode_times_and_vinf(x_best))

if __name__ == '__main__':
    main()
