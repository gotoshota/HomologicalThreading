import numpy as np
import matplotlib.pyplot as plt
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", nargs="+", type=str, help="filename")
    return parser.parse_args()


def main():
    args = get_args()
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 1)
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot(axes[i], filename)
    plt.show()


def plot(ax, filename):
    data = np.load(filename)
    alphas = data["alphas"]
    betti = data["betti_pd_i"]
    const = 100.0
    ax.plot(alphas, np.array(betti, dtype=np.float32) / const, label="i")
    # betti = data["betti_pd_i_cup_j"]
    # ax.plot(alphas, betti / 10000, label="i cup j")
    betti = data["betti_threading"]
    const = 1.0
    print (betti[10:30] / const)
    ax.plot(alphas, betti / const, label="threading")

    # ax.set_xlim(0.1, 5000)
    ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.legend()


if __name__ == "__main__":
    main()
