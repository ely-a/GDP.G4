import pykep as pk
import pygmo as pg
import numpy as np
import pandas as pd
from datetime import datetime
from time import time
from pykep.planet import jpl_lp
from pykep.trajopt import mga


def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days


def run_optimization(algo, algo_name, generations=20, seed=None):
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2031, 1, 1)), to_mjd2000(datetime(2034, 1, 1))]
    tof = [[400, 800], [2375, 3000]]
    vinf = 11.0

    udp = mga(seq=seq, t0=t0, tof=tof, vinf=vinf, multi_objective=False, tof_encoding='direct')
    prob = pg.problem(udp)

    if seed is not None:
        algo.set_seed(seed)

    archi = pg.archipelago(algo=algo, prob=prob, n=8, pop_size=30)

    start = time()
    archi.evolve(generations)
    archi.wait()
    elapsed = time() - start

    sols_f = archi.get_champions_f()
    best_dv = min(s[0] for s in sols_f)

    return best_dv, elapsed


def main():
    algorithms = [
        (pg.algorithm(pg.sade(gen=200)), "SADE"),
        (pg.algorithm(pg.de(gen=200)), "DE"),
        (pg.algorithm(pg.pso(gen=200)), "PSO"),
        (pg.algorithm(pg.cmaes(gen=200)), "CMA-ES"),
        (pg.algorithm(pg.sga(gen=200)), "SGA"),
    ]

    print(f"{'Algorithm':<10} | {'Best ΔV_T (km/s)':>18} | {'Time to Converge (s)':>22}")
    print("-" * 56)

    for algo, name in algorithms:
        dv, t = run_optimization(algo, name)
        print(f"{name:<10} | {dv:>18.4f} | {t:>22.2f}")



if __name__ == "__main__":
    main()
