import matplotlib.pyplot as plt
import sys
import glob
import pathlib
import time
import argparse
import numpy as np

sys.path.append(str(pathlib.Path(__file__).resolve().parent.parent / "src"))
import homological_threading as ht


def get_args():
    parser = argparse.ArgumentParser(description="Homological threading")
    parser.add_argument(
        "-i", "--input", type=str, nargs="+", help="LAMMPS data file", required=True
    )
    parser.add_argument(
        "-o", "--outputdir", type=str, help="Output directory", required=True
    )
    return parser.parse_args()


def calc_ensemble_betti_numbers(pds, normalization=1.0):
    """
    Calculate ensemble Betti numbers from persistence diagrams.

    args:
    pds: list of ndarray shape (npoints, 2)
        Persistence diagrams.

    normalization: float

    returns:
    alphas: ndarray shape (npoints,)
        Alphas.

    mean_betti_numbers: ndarray shape (npoints,)
    """
    alpha, betti_number = ht.compute_betti_number(pds)
    return alpha, betti_number / normalization


def plot_betti_numbers(alphas, betti_numbers):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(alphas, betti_numbers)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\beta$")
    plt.show()


def _threading(args):
    elapsed_times = [[], [], []]  # pd_i, pd_i_cup_j, threading
    max_alpha = 10000
    delta_alpha = 0.2
    for filename in args.input:
        # /path/to/xxx.data -> xxx.h5
        outputFile = pathlib.Path(filename).stem + ".h5"
        output_path = pathlib.Path(args.outputdir) / outputFile

        # Single chain
        pds = ht.HomologicalThreading()
        coords = pds.read_lmpdata(filename)
        time_start = time.time()
        pds.pd_i.compute(coords, dim=1, mp=False)
        betti = pds.pd_i.betti(max_alpha, delta_alpha)
        time_end = time.time()
        elapsed_times[0].append(time_end - time_start)

        # Pair of chains
        time_start = time.time()
        pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
        betti = pds.pd_i_cup_j.betti(max_alpha, delta_alpha)
        time_end = time.time()
        elapsed_times[1].append(time_end - time_start)

        # Threading
        time_start = time.time()
        pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
        betti = pds.threading.betti(max_alpha, delta_alpha)
        time_end = time.time()
        elapsed_times[2].append(time_end - time_start)

        # Save persistence diagrams
        pds.to_hdf5(output_path)

    print("Mean elapsed time for computing pd_i: ", np.mean(elapsed_times[0]))
    print("Mean elapsed time for computing pd_i_cup_j: ", np.mean(elapsed_times[1]))
    print("Mean elapsed time for computing threading: ", np.mean(elapsed_times[2]))


def _betti(args):
    outputFile = "betti.h5"
    output_path = pathlib.Path(args.outputdir) / outputFile
    pds = ht.HomologicalThreading()
    max_alpha = 5000
    delta_alpha = 0.1
    betti_pd_i = []
    betti_pd_i_cup_j = []
    betti_threading = []
    for filename in args.input:
        pds.from_hdf5(filename)
        alphas, betti = pds.pd_i.betti(max_alpha, delta_alpha)
        betti_pd_i.append(betti)
        alphas, betti = pds.pd_i_cup_j.betti(max_alpha, delta_alpha)
        betti_pd_i_cup_j.append(betti)
        alphas, betti = pds.threading.betti(max_alpha, delta_alpha)
        betti_threading.append(betti)
    betti_pd_i = np.array(betti_pd_i)
    betti_pd_i_cup_j = np.array(betti_pd_i_cup_j)
    betti_threading = np.array(betti_threading)
    alphas = np.array(alphas)

    betti_pd_i_mean = np.mean(betti_pd_i, axis=0)
    betti_pd_i_cup_j_mean = np.mean(betti_pd_i_cup_j, axis=0)
    betti_threading_mean = np.mean(betti_threading, axis=0)

    np.savez(
        output_path,
        alphas=alphas,
        betti_pd_i=betti_pd_i_mean,
        betti_pd_i_cup_j=betti_pd_i_cup_j_mean,
        betti_threading=betti_threading_mean,
    )


def main():
    args = get_args()
    # _threading(args)
    _betti(args)


if __name__ == "__main__":
    main()
