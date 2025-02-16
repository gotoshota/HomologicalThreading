import matplotlib.pyplot as plt
import sys
import pathlib
import time

sys.path.append(str(pathlib.Path(__file__).resolve().parent.parent / "src"))

# from homological_threading import HomologicalThreading
import homological_threading as ht


def plot_init():
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "serif"

    figs, ax = plt.subplots()
    return figs, ax


def plot_betti_numbers(pd, normalization=1.0, title=None):
    """ 
    Plot Betti numbers from persistence diagrams.

    args:
    pd: ndarray shape (npoints, 2)
        Persistence diagrams.

    normalization: float
        Normalization factor. Default is 1.0.

    title: str
    """
    betti_numbers = []
    alphas, betti_number = ht.compute_betti_number(pd)

    fig, ax = plot_init()
    ax.plot(alphas, betti_number / normalization)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\langle \beta_1 \rangle$")
    ax.set_xscale("log")
    if title is not None:
        ax.set_title(title)
    plt.show()


def test(filename, output):
    time_start = time.time()
    pds = ht.HomologicalThreading()
    coords = pds.read_lmpdata(filename)
    pds.pd_i.compute(coords, dim=1, mp=False)
    time_end = time.time()
    print("Elapsed time for computing pd_i: ", time_end - time_start)
    time_start = time.time()
    pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
    time_end = time.time()
    print("Elapsed time for computing pd_i_cup_j: ", time_end - time_start)
    time_start = time.time()
    pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
    time_end = time.time()
    print("Elapsed time for computing threading: ", time_end - time_start)
    pds.to_hdf5(output)
    # Test Betti numbers
    nchains = pds.pd_i.pd.shape[0]
    all_pds = pds.pd_i.pd
    # reshape: (nchains, npoints, 2) -> (nchains * npoints, 2)
    nchains = all_pds.shape[0]
    all_pds = all_pds.reshape(-1, 2)
    plot_betti_numbers(
        all_pds,
        normalization=nchains,
        title=r"ensemble averaged Betti number from PD$(i)$",
    )
    all_pds = pds.pd_i_cup_j.pd
    # rehsape: (nchains, nchains, npoints, 2) -> (nchains * nchains * npoints, 2)
    all_pds = all_pds.reshape(-1, 2)
    plot_betti_numbers(
        all_pds,
        normalization=nchains * (nchains - 1),
        title=r"ensemble averaged Betti number from PD$(i) \cup$ PD$(j)$",
    )
    all_pds = pds.threading.pd
    # reshape: (nchains, nchains, npoints, 2) -> (nchains * nchains * npoints, 2)
    all_pds = all_pds.reshape(-1, 2)
    plot_betti_numbers(
        all_pds,
        normalization=nchains * (nchains - 1),
        title=r"ensemble averaged Betti number from threading",
    )


if __name__ == "__main__":
    test("../data/N10M100.data", "../data/ht.h5")
