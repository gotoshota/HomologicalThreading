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
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # PD command
    pd_parser = subparsers.add_parser("pd", help="Compute persistence diagrams")
    pd_parser.add_argument("-i", "--input", nargs="+", help="Input LAMMPS DATA files")
    pd_parser.add_argument("-o", "--outputdir", default=".", help="Output directory")

    # Betti command
    betti_parser = subparsers.add_parser("betti", help="Compute Betti numbers")
    betti_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files")
    betti_parser.add_argument("-o", "--outputdir", default=".", help="Output directory")

    # Num threading command
    num_threading_parser = subparsers.add_parser("num_threading", help="Number of threading")
    num_threading_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files")

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
        time_end = time.time()
        elapsed_times[0].append(time_end - time_start)

        # Pair of chains
        time_start = time.time()
        pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
        time_end = time.time()
        elapsed_times[1].append(time_end - time_start)

        # Threading
        time_start = time.time()
        pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
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
    betti_pd_i = np.zeros(int(max_alpha / delta_alpha + 1))
    betti_pd_i_cup_j = np.zeros(int(max_alpha / delta_alpha + 1))
    betti_threading = np.zeros(int(max_alpha / delta_alpha + 1))
    for filename in args.input:
        pds.from_hdf5(filename)
        alphas, betti = pds.pd_i.betti(max_alpha, delta_alpha)
        betti_pd_i += betti
        alphas, betti = pds.pd_i_cup_j.betti(max_alpha, delta_alpha)
        betti_pd_i_cup_j += betti
        alphas, betti = pds.threading.betti(max_alpha, delta_alpha)
        betti_threading += betti

    betti_pd_i /= len(args.input)
    betti_pd_i_cup_j /= len(args.input)
    betti_threading /= len(args.input)
    np.savez(
        output_path,
        alphas=alphas,
        betti_pd_i=betti_pd_i,
        betti_pd_i_cup_j=betti_pd_i_cup_j,
        betti_threading=betti_threading,
    )

def _num_threading(args):
    for filename in args.input:
        pds = ht.HomologicalThreading()
        pds.from_hdf5(filename)
        n_a, n_p = pds.threading.num_threading()
        print(pds.threading.flags[0])
        print(n_a)
        print(np.mean(n_a))
        print(np.std(n_a))
        print(n_p)
        print(np.mean(n_p))
        print(np.std(n_p))

def main():
    args = get_args()
    if args.command == "pd":
        _threading(args)
    elif args.command == "betti":
        _betti(args)
    elif args.command == "num_threading":
        _num_threading(args)


if __name__ == "__main__":
    main()
