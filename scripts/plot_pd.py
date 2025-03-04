import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", nargs="+", type=str, help="filename")
    return parser.parse_args()


def main():
    args = get_args()
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 3, figsize=(15, 5*num_files))
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot(axes[i], filename)
    # All axes should be [0, 400]
    for ax in axes:
        for a in ax:
            a.set_xlim([0, 50])
            a.set_ylim([0, 50])
    plt.savefig("pd.png")


def plot(ax, filename):
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




if __name__ == "__main__":
    main()

