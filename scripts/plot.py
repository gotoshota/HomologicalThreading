import numpy as np
import matplotlib.pyplot as plt
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # PD command
    pd_parser = subparsers.add_parser("pd", help="Compute persistence diagrams")
    pd_parser.add_argument("-i", "--input", nargs="+", help="Input LAMMPS DATA files")
    pd_parser.add_argument("-o", "--output", default=None, help="Output image file")
    
    # Betti command
    betti_parser = subparsers.add_parser("betti", help="Compute Betti numbers")
    betti_parser.add_argument("-i", "--input", nargs="+", help="Input npz files")
    betti_parser.add_argument("-o", "--output", default=None, help="Output image file")

    return parser.parse_args()


def betti(args):
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 1)
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot(axes[i], filename)
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()


def plot_betti(ax, filename):
    data = np.load(filename)
    alphas = data["alphas"]
    betti = data["betti_pd_i"]
    const = 100.0
    ax.plot(np.sqrt(alphas), np.array(betti, dtype=np.float32) / const, label="i")
    # betti = data["betti_pd_i_cup_j"]
    # ax.plot(alphas, betti / 10000, label="i cup j")
    betti = data["betti_threading"]
    const = 1.0
    print (betti[10:30] / const)
    ax.plot(np.sqrt(alphas), betti / const, label="threading")

    ax.set_xlim(0, 10)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.legend()

def pd(args):
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 3, figsize=(15, 5*num_files))
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot_pd(axes[i], filename)
    # All axes should be [0, 400]
    for ax in axes:
        for a in ax:
            a.set_xlim([0, 50])
            a.set_ylim([0, 50])
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()

def plot_pd(ax, filename):
    with h5py.File(filename, "r") as data:
        tmp = data["pd_i/pd"][:]
        tmp = tmp.reshape(-1, 2)
        # #. of points without nan
        print(np.sum(~np.isnan(tmp[:, 0])))
        ax[0].scatter(tmp[:, 0], tmp[:, 1])
        ax[0].set_title("pd_i")
        tmp = data["pd_i_cup_j/pd"][:]
        tmp = tmp.reshape(-1, 2)
        print(np.sum(~np.isnan(tmp[:, 0])))
        ax[1].scatter(tmp[:, 0], tmp[:, 1])
        ax[1].set_title("pd_i_cup_j")
        tmp = data["threading/pd"][:]
        tmp = tmp.reshape(-1, 2)
        print(np.sum(~np.isnan(tmp[:, 0])))
        ax[2].scatter(tmp[:, 0], tmp[:, 1])
        ax[2].set_title("threading")


def main(args):
    if args.command == "pd":
        pd(args)
    elif args.command == "betti":
        betti(args)
    else:
        raise ValueError("Unknown command")



if __name__ == "__main__":
    args = get_args()
    main()
