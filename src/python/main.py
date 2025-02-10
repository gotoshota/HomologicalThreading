import lammps_io as io  # my code
import homcloud.interface as hc
import numpy as np
import h5py


class HomologicalThreading:
    """
    Class for computing the homological threading of ring polymers.
    """

    def __init__(self):
        self.pd_i = self.PD_i()
        self.pd_i_cup_j = self.PD_i_cup_j()
        self.pd_threading = self.PD_threading()

    class PD_i:
        """
        Class for storing the persistence diagram of single ring polymer.
        """

        def __init__(self):
            # shape: (nchains, npoints, 2) 0: birth, 1: death
            self.pd = None  # 最終的には np.array に変換する

        def compute(self, coords, dim=1):
            """
            Compute the persistence diagram of a single ring polymer.

            args:
            coords: np.array, shape=(nchains, nbeads, 3)
            nbeads: int, number of beads in the polymer
            nchains: int, number of chains in the polymer
            dim: int, dimension of the homology group to compute
            """
            nchains = coords.shape[0]
            pd_list = []  # 各チェインのPDを格納するリスト (shape: (npoints, 2))
            for i in range(nchains):
                polymer_coords = coords[i]
                tmp = hc.PDList.from_alpha_filtration(polymer_coords)
                pd_obj = tmp.dth_diagram(dim)
                pd_chain = np.array(
                    [pd_obj.births, pd_obj.deaths]
                ).T  # (npoints, 2) 0: birth, 1: death
                pd_list.append(pd_chain)
            # pd_list [shape: (nchains, npoints, 2)] の npoints を揃えるため，padding 処理を行う
            max_npoints = max([len(pd) for pd in pd_list])
            # padding した空の配列を作成
            pd_array = np.full((nchains, max_npoints, 2), np.nan)
            # pd_list の各要素を pd_array にコピー
            for i, pd in enumerate(pd_list):
                pd_array[i, : len(pd)] = np.array(pd)
            self.pd = pd_array

    class PD_i_cup_j:
        """
        Class for storing the persistence diagram of the cup product of two ring polymers.
        """

        def __init__(self):
            # shape: (nchains, nchains, npoints, 2) 0: birth, 1: death
            self.pd = None

        def compute(self, coords, dim=1):
            """
            Compute the persistence diagram of the cup product of two ring polymers.

            args:
            coords: np.array, shape=(nchains, nbeads, 3)
            nbeads: int, number of beads in the polymer
            nchains: int, number of chains in the polymer
            dim: int, dimension of the homology group to compute
            """
            nchains = coords.shape[0]
            pd_list = []  # shape: (nchains, nchains, npoints, 2)
            for i in range(nchains):
                pd_list_chain = []  # shape: (nchains, npoints, 2)
                for j in range(nchains):
                    if i == j:
                        pd_list_chain.append([np.nan, np.nan])
                        continue  # 同じチェイン同士は計算しない
                    polymer_coords_i = coords[i]
                    polymer_coords_j = coords[j]
                    # Compute the persistence diagram of the cup product
                    tmp = hc.PDList.from_alpha_filtration(
                        np.concatenate([polymer_coords_i, polymer_coords_j])
                    )
                    pd_obj = tmp.dth_diagram(dim)
                    pd_chain = np.array(
                        [pd_obj.births, pd_obj.deaths]
                    ).T  # shape: (npoints, 2)
                    pd_list_chain.append(pd_chain)  # shape: (nchains, npoints, 2)
                pd_list.append(pd_list_chain)
            # pd_list [shape: (nchains, nchains, npoints, 2)] の npoints を揃えるため，padding 処理を行う
            max_npoints = max(
                [len(pd) for pd_list_chain in pd_list for pd in pd_list_chain]
            )
            # padding した空の配列を作成
            pd_array = np.full((nchains, nchains, max_npoints, 2), np.nan)
            # pd_list の各要素を pd_array にコピー
            for i, pd_list_chain in enumerate(pd_list):
                for j, pd in enumerate(pd_list_chain):
                    pd_array[i, j, : len(pd)] = np.array(pd)
            self.pd = pd_array

    class PD_threading:
        """
        Class for storing the homological threading of ring polymers.
        """

        def __init__(self):
            self.pd = None

    def to_hdf5(self, filename="pd.h5"):
        """
        Save the persistence diagrams to a HDF5 file.
        """
        with h5py.File(filename, "w") as f:
            if self.pd_i.pd is not None:
                f.create_dataset("pd_i", data=self.pd_i.pd)
            if self.pd_i_cup_j.pd is not None:
                f.create_dataset("pd_i_cup_j", data=self.pd_i_cup_j.pd)
            # f.create_dataset("pd_threading", data=self.pd_threading.pd)


if __name__ == "__main__":
    import time

    # Load the coordinates of the ring polymers
    data = io.LammpsData("../../data/N10M100.data")
    coords = np.array(data.atoms.coords)
    # (nparticles, 3) -> (nchains, nbeads, 3)
    nchains = data.atoms.num_mols
    nbeads = len(coords) // nchains
    coords = coords.reshape(nchains, nbeads, 3)

    # Compute the persistence diagram of a single ring polymer
    time_start = time.time()
    pds = HomologicalThreading()
    pds.pd_i.compute(coords, dim=1)
    time_end = time.time()
    print("Elapsed time for computing pd_i: ", time_end - time_start)
    time_start = time.time()
    pds.pd_i_cup_j.compute(coords, dim=1)
    time_end = time.time()
    print("Elapsed time for computing pd_i_cup_j: ", time_end - time_start)
    pds.to_hdf5()
