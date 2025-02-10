import lammps_io as io  # my code
import homcloud.interface as hc
import numpy as np
import h5py
import multiprocessing as mp
import os


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

        def compute(self, coords, dim=1, mp=False, num_processes=None):
            """
            Compute the persistence diagram of a single ring polymer.

            args:
            coords: np.array, shape=(nchains, nbeads, 3)
            nbeads: int, number of beads in the polymer
            nchains: int, number of chains in the polymer
            dim: int, dimension of the homology group to compute
            num_processes: int, number of processes to use for parallel computation
            """
            if mp:
                self.compute_mp(coords, dim, num_processes)
            else:
                self.compute_single(coords, dim)

        def compute_single(self, coords, dim=1):
            nchains = coords.shape[0]
            pd_list = []  # 各チェインのPDを格納するリスト (shape: (npoints, 2))
            for i in range(nchains):
                polymer_coords = coords[i]
                tmp = hc.PDList.from_alpha_filtration(polymer_coords)
                pd_obj = tmp.dth_diagram(dim)
                pd_chain = np.array([pd_obj.births, pd_obj.deaths]).T
                pd_list.append(pd_chain)
            # pd_list [shape: (nchains, npoints, 2)] の npoints を揃えるため，padding 処理を行う
            max_npoints = max([len(pd) for pd in pd_list])
            # padding した空の配列を作成
            pd_array = np.full((nchains, max_npoints, 2), np.nan)
            # pd_list の各要素を pd_array にコピー
            for i, pd in enumerate(pd_list):
                pd_array[i, : len(pd)] = np.array(pd)
            self.pd = pd_array

        def compute_mp(self, coords, dim=1, num_processes=None):
            nchains = coords.shape[0]
            pd_list = []  # 各チェインのPDを格納するリスト (shape: (npoints, 2))
            # 並列プロセス数を取得
            if num_processes is None:
                num_processes = int(os.environ.get("OMP_NUM_THREADS", mp.cpu_count()))
            pool = mp.Pool(num_processes)
            # 各プロセスに割り当てる範囲を計算
            chunk_size = (nchains + num_processes - 1) // num_processes
            tasks = []
            for i in range(num_processes):
                ista = i * chunk_size
                iend = min(ista + chunk_size, nchains)
                tasks.append((ista, iend, coords, dim))
            results = pool.map(self._worker, tasks)
            pool.close()
            pool.join()
            # 並列化処理終了
            # PD のリストを結合
            pd_list.extend([pd for pd_chain in results for pd in pd_chain])
            # pd_list [shape: (nchains, npoints, 2)] の npoints を揃えるため，padding 処理を行う
            max_npoints = max([len(pd) for pd in pd_list])
            # padding した空の配列を作成
            pd_array = np.full((nchains, max_npoints, 2), np.nan)
            # pd_list の各要素を pd_array にコピー
            for i, pd in enumerate(pd_list):
                pd_array[i, : len(pd)] = np.array(pd)
            self.pd = pd_array

        def _worker(self, args):
            """
            指定された範囲の計算を行う
            """
            ista, iend, coords, dim = args
            partial_pd_list = []
            for i in range(ista, iend):
                polymer_coords = coords[i]
                tmp = hc.PDList.from_alpha_filtration(polymer_coords)
                pd_obj = tmp.dth_diagram(dim)
                pd_chain = np.array([pd_obj.births, pd_obj.deaths]).T
                partial_pd_list.append(pd_chain)
            return partial_pd_list

    class PD_i_cup_j:
        """
        Class for storing the persistence diagram of the cup product of two ring polymers.
        """

        def __init__(self):
            # shape: (nchains, nchains, npoints, 2) 0: birth, 1: death
            self.pd = None

        def compute(self, coords, dim=1, mp=False, num_processes=None):
            """
            Compute the persistence diagram of the cup product of two ring polymers.

            args:
            coords: np.array, shape=(nchains, nbeads, 3)
            nbeads: int, number of beads in the polymer
            nchains: int, number of chains in the polymer
            dim: int, dimension of the homology group to compute
            num_processes: int, number of processes to use for parallel computation
            """
            if mp:
                self.compute_mp(coords, dim, num_processes)
            else:
                self.compute_single(coords, dim)

        def compute_single(self, coords, dim=1):
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

        def compute_mp(self, coords, dim=1, num_processes=None):
            nchains = coords.shape[0]
            pd_list = []
            # 並列プロセス数を取得
            if num_processes is None:
                num_processes = int(os.environ.get("OMP_NUM_THREADS", mp.cpu_count()))
            pool = mp.Pool(num_processes)
            # 各プロセスに割り当てる範囲を計算
            chunk_size = (nchains + num_processes - 1) // num_processes
            tasks = []
            for i in range(num_processes):
                ista = i * chunk_size
                iend = min(ista + chunk_size, nchains)
                tasks.append((ista, iend, coords, dim))
            results = pool.map(self._worker, tasks)
            pool.close()
            pool.join()
            # 並列化処理終了
            # PD のリストを結合
            pd_list.extend([pd for pd_chain in results for pd in pd_chain])
            # pd_list [shape: (nchains, nchains, npoints, 2)] の npoints を揃えるため，padding 処理を行う
            max_npoints = max([len(pd) for pd_chain in pd_list for pd in pd_chain])
            # padding した空の配列を作成
            pd_array = np.full((nchains, nchains, max_npoints, 2), np.nan)
            # pd_list の各要素を pd_array にコピー
            for i, pd_list_chain in enumerate(pd_list):
                for j, pd in enumerate(pd_list_chain):
                    pd_array[i, j, : len(pd)] = np.array(pd)
            self.pd = pd_array

        def _worker(self, args):
            """
            指定された範囲の計算を行う
            """
            ista, iend, coords, dim = args
            nchains = coords.shape[0]
            partial_pd_list = []
            for i in range(ista, iend):
                pd_list_chain = []
                for j in range(nchains):
                    if i == j:
                        pd_list_chain.append([np.nan, np.nan])
                        continue
                    polymer_coords_i = coords[i]
                    polymer_coords_j = coords[j]
                    # Compute the persistence diagram of the cup product
                    tmp = hc.PDList.from_alpha_filtration(
                        np.concatenate([polymer_coords_i, polymer_coords_j])
                    )
                    pd_obj = tmp.dth_diagram(dim)
                    pd_chain = np.array([pd_obj.births, pd_obj.deaths]).T
                    pd_list_chain.append(pd_chain)
                partial_pd_list.append(pd_list_chain)
            return partial_pd_list

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
            if self.pd_threading.pd is not None:
                f.create_dataset("pd_threading", data=self.pd_threading.pd)


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
    pds.pd_i.compute(coords, dim=1, mp=True)
    time_end = time.time()
    print("Elapsed time for computing pd_i: ", time_end - time_start)
    time_start = time.time()
    pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
    time_end = time.time()
    print("Elapsed time for computing pd_i_cup_j: ", time_end - time_start)
    pds.to_hdf5()
