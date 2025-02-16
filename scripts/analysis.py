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
    parser.add_argument("-i", "--input", type=str, nargs="+", help="LAMMPS data file")
    parser.add_argument("-o", "--outputdir", type=str, help="Output directory")
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


def main():
    args = get_args()
    print(args.input)
    elapsed_times = [[], [], []]  # pd_i, pd_i_cup_j, threading
    all_pds = [[], [], []]  # pd_i, pd_i_cup_j, threading: Each shape (npoints, 2)
    for filename in args.input:
        # /path/to/xxx.data -> xxx.h5
        outputFile = pathlib.Path(filename).stem + ".h5"
        output_path = pathlib.Path(args.outputdir) / outputFile

        # debug
        pds = ht.HomologicalThreading()
        pds.from_hdf5(output_path)

        pds = ht.HomologicalThreading()
        coords = pds.read_lmpdata(filename)
        time_start = time.time()
        pds.pd_i.compute(coords, dim=1, mp=False)
        time_end = time.time()
        elapsed_times[0].append(time_end - time_start)

        time_start = time.time()
        pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
        time_end = time.time()
        elapsed_times[1].append(time_end - time_start)

        time_start = time.time()
        pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
        time_end = time.time()
        elapsed_times[2].append(time_end - time_start)

        pds.to_hdf5(output_path)

        # add pds to all_pds list: shape (nchains, npoints, 2) -> (nchains * npoints, 2)
        nchains = len(pds.pd_i.pd)
        tmp = np.reshape(pds.pd_i.pd, (-1, 2))
        all_pds[0].extend(tmp)
        tmp = np.reshape(pds.pd_i_cup_j.pd, (-1, 2))
        all_pds[1].extend(tmp)
        tmp = np.reshape(pds.threading.pd, (-1, 2))
        all_pds[2].extend(tmp)

    print("Mean elapsed time for computing pd_i: ", np.mean(elapsed_times[0]))
    print("Mean elapsed time for computing pd_i_cup_j: ", np.mean(elapsed_times[1]))
    print("Mean elapsed time for computing threading: ", np.mean(elapsed_times[2]))

    pd_i = np.array(all_pds[0])
    pd_i_cup_j = np.array(all_pds[1])
    threading = np.array(all_pds[2])
    # Calculate ensemble Betti numbers
    betti_numbers = [[], [], []]
    alphas = [[], [], []]
    normalization = nchains
    alphas[0], betti_numbers[0] = calc_ensemble_betti_numbers(pd_i, normalization)
    normalization = nchains * (nchains - 1)
    alphas[1], betti_numbers[1] = calc_ensemble_betti_numbers(pd_i_cup_j, normalization)
    normalization = nchains * (nchains - 1)
    alphas[2], betti_numbers[2] = calc_ensemble_betti_numbers(threading, normalization)

    for i in range(3):
        plot_betti_numbers(alphas[i], betti_numbers[i])


if __name__ == "__main__":
    main()
