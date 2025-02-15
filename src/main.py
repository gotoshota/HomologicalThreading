import lammps_io as io
from fortran import compute as fc
import homcloud.interface as hc
import numpy as np
import h5py
import multiprocessing as mp
import os
import sys
import time

from scipy.spatial import KDTree

if sys.stdin.closed:
    sys.stdin = open(os.devnull, "r")
if sys.stdout.closed:
    sys.stdout = open(os.devnull, "w")
if sys.stderr.closed:
    sys.stderr = open(os.devnull, "w")


version = "0.1.0"


class HomologicalThreading:
    """
    Class for computing the homological threading of ring polymers.
    """

    def __init__(self, rho=None, epsilon_theta=None, source=None):
        self.pd_i = self.PD_i(self)
        self.pd_i_cup_j = self.PD_i_cup_j(self)
        self.threading = self.Threading(self)
        self.metadata = {
            "description": "Homological threading of ring polymers",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": version,
            "nchains": None,
            "nbeads": None,
            "nparticles": None,
            "box_dim": None,
            "rho": rho,
            "epsilon_theta": epsilon_theta,
            "source": source if source else "Unknown",
            "threading_threshold": None,
        }

    def print_metadata(self):
        for key, value in self.metadata.items():
            print(f"{key}: {value}")

    def read_lmpdata(self, filename):
        """
        Read the coordinates of the ring polymers from a LAMMPS data file.

        args:
        filename: str, path to the LAMMPS data file

        return:
        coords: np.array, shape=(nchains, nbeads, 3)
        """
        data = io.LammpsData(filename)
        data.polyWrap()
        coords = np.array(data.atoms.coords)
        nchains = data.atoms.num_mols
        nbeads = data.atoms.num_atoms // nchains
        box_dim = data.box.lx

        # Metadata に情報を格納
        self.metadata["nchains"] = nchains
        self.metadata["nbeads"] = nbeads
        self.metadata["nparticles"] = nchains * nbeads
        self.metadata["box_dim"] = box_dim
        self.metadata["source"] = filename

        # (nparticles, 3) -> (nchains, nbeads, 3)
        coords = coords.reshape(nchains, nbeads, 3)
        return coords

    class PD_i:
        """
        Class for storing the persistence diagram of single ring polymer.
        """

        def __init__(self, parent):
            # shape: (nchains, npoints, 2) 0: birth, 1: death
            self.parent = parent
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
            nchains = coords.shape[0]  # coords: (nchains, nbeads, 3)
            nbeads = coords.shape[1]
            self.parent.metadata["nchains"] = nchains
            self.parent.metadata["nbeads"] = nbeads
            self.parent.metadata["nparticles"] = nchains * nbeads

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
                rows, cols = pd.shape
                pd_array[i, :rows, :cols] = pd
            self.pd = pd_array

        def compute_mp(self, coords, dim=1, num_processes=None):
            nchains = coords.shape[0]
            nbeads = coords.shape[1]
            self.parent.metadata["nchains"] = nchains
            self.parent.metadata["nbeads"] = nbeads
            self.parent.metadata["nparticles"] = nchains * nbeads
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

        def __init__(self, parent):
            # shape: (nchains, nchains, npoints, 2) 0: birth, 1: death
            self.parent = parent
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
            nbeads = coords.shape[1]
            self.parent.metadata["nchains"] = nchains
            self.parent.metadata["nbeads"] = nbeads
            self.parent.metadata["nparticles"] = nchains * nbeads
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
            nbeads = coords.shape[1]
            self.parent.metadata["nchains"] = nchains
            self.parent.metadata["nbeads"] = nbeads
            self.parent.metadata["nparticles"] = nchains * nbeads
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

        def __init__(self, parent):
            self.parent = parent
            self.flags = None  # shape: (nchains, nchains, npoints same as pd_i)
            self.pd = None  # shape: (nchains, nchains, npoints, 2)

        def compute(self, pd_i, pd_i_cup_j, threshold=1e-10):
            """
            Compute the homological threading of ring polymers.

            args:
            pd_i: np.array, shape=(nchains, npoints, 2)
            pd_i_cup_j: np.array, shape=(nchains, nchains, npoints', 2)
            """
            nchains = pd_i.shape[0]
            npoints = pd_i.shape[1]
            self.parent.metadata["threading_threshold"] = threshold

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
            fc.threading(pd_i_fort, pd_i_cup_j_fort, flags_fort, pd_fort, threshold)
            mask = pd_fort == -1
            pd_fort[mask] = np.nan

            # Fortran からの返り値を python 用に変換
            self.pd = pd_fort.T
            self.flags = flags_fort.astype(bool)

        def compute_kdtree(self, pd_i, pd_i_cup_j, tol=1e-10):
            """
            Copute the homological threading of ring polymers using KDTree.

            args:
            pd_i: np.array, shape=(nchains, npoints, 2)
            pd_i_cup_j: np.array, shape=(nchains, nchains, npoints', 2)
            tol: float, tolerance for the KDTree
            """
            print("This function is not implemented yet.")
            print("Please use the compute method.")
            return

            nchains = pd_i.shape[0]
            npoints = pd_i.shape[1]
            threading_pd = np.zeros((nchains, nchains, npoints, 2), dtype=np.float64)
            for passive_idx in range(nchains):
                for active_idx in range(nchains):
                    threading_pd[passive_idx, active_idx] = self.set_difference(
                        pd_i[passive_idx], pd_i_cup_j[passive_idx, active_idx], tol
                    )

        def set_difference(self, A, B, tol=1e-10):
            """
            Compute the difference of two persistence diagrams.

            args:
            A: np.array, shape=(npoints, 2)
            B: np.array, shape=(npoints', 2)
            tol: float, tolerance for the KDTree

            return C = A / B
            """
            if B.size == 0:
                return A.copy()

            # Bに対して KDTree を構築
            # nan は 削除する
            B_prime = B[~np.isnan(B).any(axis=1)]
            if B_prime.size == 0:
                return A.copy()
            tree = KDTree(B_prime)
            # Aの各点に対して最近傍点を探索
            C_prime = []
            for point in A:
                if np.isnan(point).any():
                    continue
                if not tree.query_ball_point(point, tol):
                    C_prime.append(point)

            if len(C_prime) == 0:
                C = np.full(A.shape, np.nan)
            else:
                # C を padding して返す
                C = np.full(A.shape, np.nan)
                C[: len(C_prime), :] = np.array(C_prime)

            return C

    def to_hdf5(self, filename):
        """
        Save the persistence diagrams to a HDF5 file.
        """
        with h5py.File(filename, "w") as f:
            if self.pd_i.pd is not None:
                f.create_group("pd_i")
                f.create_dataset("pd_i/pd", data=self.pd_i.pd)
            if self.pd_i_cup_j.pd is not None:
                f.create_group("pd_i_cup_j")
                f.create_dataset("pd_i_cup_j/pd", data=self.pd_i_cup_j.pd)
            if self.threading.flags is not None or self.threading.pd is not None:
                f.create_group("threading")
                if self.threading.flags is not None:
                    f.create_dataset("threading/flags", data=self.threading.flags)
                if self.threading.pd is not None:
                    f.create_dataset("threading/pd", data=self.threading.pd)
            f.create_group("Metadata")
            for key, value in self.metadata.items():
                if value is None:
                    value = "None"
                f["Metadata"].attrs[key] = value

    def from_hdf5(self, filename):
        """
        Load the persistence diagrams from a HDF5 file.
        """
        with h5py.File(filename, "r") as f:
            if "pd_i" in f:
                self.pd_i.pd = f["pd_i/pd"][:]
            if "pd_i_cup_j" in f:
                self.pd_i_cup_j.pd = f["pd_i_cup_j/pd"][:]
            if "threading" in f:
                self.threading.flags = f["threading/flags"][:]
                self.threading.pd = f["threading/pd"][:]
            if "Metadata" in f:
                meta_grp = f["Metadata"]
                for key in meta_grp.attrs:
                    self.metadata[key] = meta_grp.attrs[key]
                    if self.metadata[key] == "None":
                        self.metadata[key] = None


def compute_betti_number(pd, d_alpha=0.01):
    """
    Compute the Betti number from the persistence diagram.

    args:
    pd: np.array, shape=(npoints, 2)

    return:
    betti_numbers: np.array, shape=(n_alpha)
    """
    # fortran 用に配列を変換
    pd_fort = pd.T
    pd_fort = np.asfortranarray(pd_fort)

    # nan を除去
    tmp = pd[~np.isnan(pd).any(axis=1)]
    max_death = np.max(tmp[:, 1])
    print(max_death)
    n_alpha = int(max_death / d_alpha) + 1
    betti_number = fc.betti_number(pd_fort, d_alpha, n_alpha)
    alphas = np.arange(0, n_alpha * d_alpha, d_alpha)
    return alphas, betti_number


def test(filename, output):
    time_start = time.time()
    pds = HomologicalThreading()
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
    betti_numbers = []
    all_pds = pds.pd_i.pd
    # reshape: (nchains, npoints, 2) -> (nchains * npoints, 2)
    nchains = all_pds.shape[0]
    all_pds = all_pds.reshape(-1, 2)
    alphas, betti_number = compute_betti_number(all_pds)
    import matplotlib.pyplot as plt
    figs, ax = plt.subplots()
    ax.plot(alphas, betti_number/nchains)
    plt.show()


if __name__ == "__main__":
    import time

    filename = sys.argv[1]
    output = sys.argv[2]
    time_start = time.time()
    test(filename, output)
    time_end = time.time()
    print("Total elapsed time: ", time_end - time_start)
