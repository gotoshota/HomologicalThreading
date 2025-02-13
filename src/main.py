import lammps_io as io  # my code
import homcloud.interface as hc
import numpy as np
import h5py
import multiprocessing as mp
import fortran.compute as fc
import os
import sys

if sys.stdin.closed:
    sys.stdin = open(os.devnull, "r")
if sys.stdout.closed:
    sys.stdout = open(os.devnull, "w")
if sys.stderr.closed:
    sys.stderr = open(os.devnull, "w")


class HomologicalThreading:
    """
    Class for computing the homological threading of ring polymers.
    """

    def __init__(self):
        self.pd_i = self.PD_i()
        self.pd_i_cup_j = self.PD_i_cup_j()
        self.threading = self.Threading()

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
            # 結果は各プロセス毎のリストになっているので、平坦化する
            all_results = [item for sublist in results for item in sublist]
            # チェインインデックスでソートする
            all_results.sort(key=lambda x: x[0])
            # ソート後、チェインごとの pd_chain 部分のみ抽出
            pd_list = [pd_chain for (_, pd_chain) in all_results]
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
            戻り値は (chain_index, pd_chain) のリスト
            """
            ista, iend, coords, dim = args
            partial_pd_list = []
            for i in range(ista, iend):
                polymer_coords = coords[i]
                tmp = hc.PDList.from_alpha_filtration(polymer_coords)
                pd_obj = tmp.dth_diagram(dim)
                pd_chain = np.array([pd_obj.births, pd_obj.deaths]).T
                partial_pd_list.append((i, pd_chain))
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
            # 結果は各プロセス毎のリストになっているので，平坦化する
            all_results = [item for sublist in results for item in sublist]
            # 外側（チェイン i）のインデックスでソート
            all_results.sort(key=lambda x: x[0])
            # 各チェイン i について，内側の結果 pd_list_chain を取り出す
            # ※ 内側のリストは _worker 内で既に j の昇順にソートしているが，念のためソートしておく
            pd_list = []
            for i, pd_list_chain in all_results:
                pd_list_chain.sort(key=lambda x: x[0])
                # 内側は (j, pd_chain) のタプルなので，pd_chain 部分のみ抽出
                pd_list.append([pd for (j, pd) in pd_list_chain])
            # pd_list の形は (nchains, nchains, variable npoints, 2)
            # 各結果の npoints（点の数）が異なるため，最大値を取得してパディングする
            max_npoints = 0
            for chain in pd_list:
                for pd in chain:
                    max_npoints = max(max_npoints, len(pd))
            # 全体の pd_array を作成（不足部分は NaN で埋める）
            pd_array = np.full((nchains, nchains, max_npoints, 2), np.nan)
            for i, chain in enumerate(pd_list):
                for j, pd in enumerate(chain):
                    pd = np.array(pd)
                    pd_array[i, j, : len(pd)] = pd
            self.pd = pd_array

        def _worker(self, args):
            """
            指定された範囲の計算を行う
            戻り値は (chain_index, pd_chain) のリスト
            """
            ista, iend, coords, dim = args
            nchains = coords.shape[0]
            partial_pd_list = []
            for i in range(ista, iend):
                pd_list_chain = []
                for j in range(nchains):
                    if i == j:
                        pd_list_chain.append((j, [np.nan, np.nan]))
                        continue
                    polymer_coords_i = coords[i]
                    polymer_coords_j = coords[j]
                    # Compute the persistence diagram of the cup product
                    tmp = hc.PDList.from_alpha_filtration(
                        np.concatenate([polymer_coords_i, polymer_coords_j])
                    )
                    pd_obj = tmp.dth_diagram(dim)
                    pd_chain = np.array([pd_obj.births, pd_obj.deaths]).T
                    pd_list_chain.append((j, pd_chain))
                # 内側のリストをチェイン j のインデックスでソート
                pd_list_chain.sort(key=lambda x: x[0])
                partial_pd_list.append((i, pd_list_chain))
            return partial_pd_list

    class Threading:
        """
        Class for storing the homological threading of ring polymers.
        """

        def __init__(self):
            self.flags = None  # shape: (nchains, nchains, npoints same as pd_i)
            self.pd = None  # shape: (nchains, nchains, npoints, 2)

        def compute(self, pd_i, pd_i_cup_j):
            """
            Compute the homological threading of ring polymers.

            args:
            pd_i: np.array, shape=(nchains, npoints, 2)
            pd_i_cup_j: np.array, shape=(nchains, nchains, npoints, 2)
            """
            nchains = pd_i.shape[0]
            npoints = pd_i.shape[1]
            # Fortran 用に配列を用意
            # pd_fort : (2, npoints, active, passive)
            pd_fort = np.zeros((2, npoints, nchains, nchains), dtype=np.float64)
            pd_fort = np.asfortranarray(pd_fort)
            # pd_i: (passive, npoints, 2) -> pd_i_fort: (2, npoints, passive)
            pd_i_fort = np.asfortranarray(pd_i.T)
            # pd_i_cup_j: (passive, active, npoints, 2) ->
            # pd_i_cup_j_fort: (2, npoints, active, passive)
            pd_i_cup_j_fort = np.asfortranarray(pd_i_cup_j.T)
            # flags_fort: (active, passive)
            flags_fort = np.ones((nchains, nchains), dtype=np.int32)
            flags_fort = np.asfortranarray(flags_fort)
            # Fortran で homological threading を計算
            fc.threading(pd_i_fort, pd_i_cup_j_fort, flags_fort, pd_fort)
            # Fortran からの返り値を python 用に変換
            pd_i = pd_i_fort.T
            pd_i_cup_j = pd_i_cup_j_fort.T
            self.pd = pd_fort.T
            self.flags = flags_fort.astype(bool)

    def to_hdf5(self, filename="pd.h5"):
        """
        Save the persistence diagrams to a HDF5 file.
        """
        with h5py.File(filename, "w") as f:
            if self.pd_i.pd is not None:
                f.create_dataset("pd_i", data=self.pd_i.pd)
            if self.pd_i_cup_j.pd is not None:
                f.create_dataset("pd_i_cup_j", data=self.pd_i_cup_j.pd)
            if self.threading.flags is not None:
                f.create_dataset("threading_flags", data=self.threading.flags)
            if self.threading.pd is not None:
                f.create_dataset("threading_pd", data=self.threading.pd)


if __name__ == "__main__":
    import time

    # Load the coordinates of the ring polymers
    # data = io.LammpsData("../data/N10M100.data")
    data = io.LammpsData("../../../Downloads/N100M100.data")
    coords = np.array(data.atoms.coords)
    # (nparticles, 3) -> (nchains, nbeads, 3)
    nchains = data.atoms.num_mols
    nbeads = len(coords) // nchains
    coords = coords.reshape(nchains, nbeads, 3)

    # Compute the persistence diagram of a single ring polymer
    print("mp=False")
    time_start = time.time()
    pds = HomologicalThreading()
    pds.pd_i.compute(coords, dim=1, mp=False)
    time_end = time.time()
    print("Elapsed time for computing pd_i: ", time_end - time_start)
    time_start = time.time()
    pds.pd_i_cup_j.compute(coords, dim=1, mp=False)
    time_end = time.time()
    print("Elapsed time for computing pd_i_cup_j: ", time_end - time_start)
    time_start = time.time()
    pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
    time_end = time.time()
    print("Elapsed time for computing threading: ", time_end - time_start)
    pds.to_hdf5("pd.h5")

    print("\nmp=True")
    time_start = time.time()
    pds = HomologicalThreading()
    pds.pd_i.compute(coords, dim=1, mp=True)
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
    pds.to_hdf5("pd_mp.h5")
