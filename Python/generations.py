import pykep as pk
import pygmo as pg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from time import time
from pykep.planet import jpl_lp
from pykep.trajopt import mga


def to_mjd2000(date):
    return (date - datetime(2000, 1, 1)).days


def run_sade_optimization(generations, seed=None):
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('neptune')]
    t0 = [to_mjd2000(datetime(2031, 1, 1)), to_mjd2000(datetime(2034, 1, 1))]
    tof = [[400, 800], [2375, 3000]]
    vinf = 11.0

    udp = mga(seq=seq, t0=t0, tof=tof, vinf=vinf, multi_objective=False, tof_encoding='direct')
    prob = pg.problem(udp)
    algo = pg.algorithm(pg.sade(gen=generations))

    if seed is not None:
        algo.set_seed(seed)

    archi = pg.archipelago(algo=algo, prob=prob, n=8, pop_size=30)

    start = time()
    archi.evolve(1)
    archi.wait()
    elapsed = time() - start

    sols_f = archi.get_champions_f()
    best_dv = min(s[0] for s in sols_f)

    return best_dv/1000, elapsed


def main():
    # Sensitivity test parameters
    gen_start = 25
    gen_end = 300
    gen_step = 25

    generation_values = list(range(gen_start, gen_end + 1, gen_step))
    results_dv = []
    results_time = []

    for gen in generation_values:
        dv, t = run_sade_optimization(gen)
        results_dv.append(dv)
        results_time.append(t)

    # Plot results
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(generation_values, results_dv, marker='o')
    plt.xlabel("Generations", fontsize=18)
    plt.ylabel("Best ΔV_T (km/s)", fontsize=18)
    plt.tick_params(axis='both', labelsize=17)
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(generation_values, results_time, marker='o', color='orange')
    plt.xlabel("Generations", fontsize=18)
    plt.ylabel("Time (s)", fontsize=18)
    plt.tick_params(axis='both', labelsize=17)
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Print differences from the final generation
    best_final_dv = results_dv[-1]
    print("\nΔV_T difference compared to highest generation result:")
    for gen, dv in zip(generation_values[:-1], results_dv[:-1]):
        abs_diff = dv - best_final_dv
        pct_diff = 100 * abs_diff / best_final_dv
        print(f"Gen {gen:>3d}: ΔV_T = {dv:.4f} km/s  |  Δ = {abs_diff:.4f} km/s ({pct_diff:.15f}%)")


if __name__ == "__main__":
    main()
