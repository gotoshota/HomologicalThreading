import matplotlib.pyplot as plt
import sys
import glob
import pathlib
import time
import argparse
import numpy as np
import h5py

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
    betti_parser.add_argument("-f", "--output_file", default="analysis.h5", help="Output HDF5 file name")

    # Num threading command
    num_threading_parser = subparsers.add_parser("num_threading", help="Number of threading")
    num_threading_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files")
    num_threading_parser.add_argument("-o", "--outputdir", default=".", help="Output directory")
    num_threading_parser.add_argument("-f", "--output_file", default="analysis.h5", help="Output HDF5 file name")

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
    output_file = args.output_file
    output_path = pathlib.Path(args.outputdir) / output_file
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
    
    # HDF5形式で保存
    # ファイルが存在する場合は追記モードで開く
    file_mode = 'a' if output_path.exists() else 'w'
    with h5py.File(str(output_path), file_mode) as f:
        # bettiグループが存在しない場合は作成、存在する場合は削除して再作成
        if 'betti' in f:
            del f['betti']
        betti_group = f.create_group('betti')
            
        # データセットを作成
        betti_group.create_dataset('alphas', data=alphas)
        betti_group.create_dataset('pd_i', data=betti_pd_i)
        betti_group.create_dataset('pd_i_cup_j', data=betti_pd_i_cup_j)
        betti_group.create_dataset('threading', data=betti_threading)
        
        # メタデータの追加
        betti_group.attrs['num_files'] = len(args.input)
        betti_group.attrs['max_alpha'] = max_alpha
        betti_group.attrs['delta_alpha'] = delta_alpha
        betti_group.attrs['timestamp'] = time.time()
        betti_group.attrs['input_files'] = ','.join([str(f) for f in args.input])

def _num_threading(args):
    output_file = args.output_file
    output_path = pathlib.Path(args.outputdir) / output_file
    nchains = None
    
    # ファイルが存在する場合は追記モードで開く
    file_mode = 'a' if output_path.exists() else 'w'
    with h5py.File(str(output_path), file_mode) as f:
        # num_threadingグループが存在しない場合は作成、存在する場合は削除して再作成
        if 'num_threading' in f:
            del f['num_threading']
        threading_group = f.create_group('num_threading')
        
        # n_aとn_pのサブグループを作成
        n_a_group = threading_group.create_group('n_a')
        n_p_group = threading_group.create_group('n_p')
        
        # 入力ファイルのメタデータを保存
        threading_group.attrs['num_files'] = len(args.input)
        threading_group.attrs['timestamp'] = time.time()
        threading_group.attrs['input_files'] = ','.join([str(f) for f in args.input])
        
        # まず最初のファイルを処理してnchainsを取得
        for filename in args.input:
            pds = ht.HomologicalThreading()
            pds.from_hdf5(filename)
            if nchains is None:
                nchains = pds.metadata['nchains']
                break
        
        # n_aとn_pのビンを設定 (0, 1, 2, ..., nchains)
        bins_n = np.arange(nchains + 2)  # 0からnchainsまでのビン + 境界用に1つ追加
        
        # 全体のヒストグラムを初期化
        hist_n_a_all = np.zeros(nchains + 1)
        hist_n_p_all = np.zeros(nchains + 1)
        
        # 全体のカウント
        total_count_n_a = 0
        total_count_n_p = 0
        
        # 全体の平均と標準偏差の計算用
        sum_n_a_mean = 0
        sum_n_a_std = 0
        sum_n_p_mean = 0
        sum_n_p_std = 0
        
        # 各ファイルのデータを処理して統計量を集計
        for filename in args.input:
            pds = ht.HomologicalThreading()
            pds.from_hdf5(filename)
            n_a, n_p = pds.threading.num_threading()
            
            # 基本統計量を計算
            n_a_mean = np.mean(n_a)
            n_a_std = np.std(n_a)
            n_p_mean = np.mean(n_p)
            n_p_std = np.std(n_p)
            
            # 全体の統計に加算
            sum_n_a_mean += n_a_mean
            sum_n_a_std += n_a_std
            sum_n_p_mean += n_p_mean
            sum_n_p_std += n_p_std
            
            # ヒストグラムを計算して全体に加算
            hist_n_a, _ = np.histogram(n_a, bins=bins_n)
            hist_n_p, _ = np.histogram(n_p, bins=bins_n)
            hist_n_a_all += hist_n_a
            hist_n_p_all += hist_n_p
            
            # カウントを加算
            total_count_n_a += len(n_a)
            total_count_n_p += len(n_p)
            
            # デバッグ出力
            print(f"処理中: {filename}")
            print(f"n_a 平均: {n_a_mean}, 標準偏差: {n_a_std}")
            print(f"n_p 平均: {n_p_mean}, 標準偏差: {n_p_std}")
        
        # 平均の平均を計算
        avg_n_a_mean = sum_n_a_mean / len(args.input) if args.input else 0
        avg_n_a_std = sum_n_a_std / len(args.input) if args.input else 0
        avg_n_p_mean = sum_n_p_mean / len(args.input) if args.input else 0
        avg_n_p_std = sum_n_p_std / len(args.input) if args.input else 0
        
        # 全体の統計情報を保存
        if total_count_n_a > 0 and total_count_n_p > 0:
            # ビン値（各ビンの中心値）を設定
            bin_centers = np.arange(nchains + 1)  # 0, 1, 2, ..., nchains
            
            # n_aの統計情報を計算して保存
            if hist_n_a_all.sum() > 0:
                pdf_n_a_all = hist_n_a_all / hist_n_a_all.sum()
                n_a_mean_all = np.sum(bin_centers * pdf_n_a_all)
                n_a_var_all = np.sum(((bin_centers - n_a_mean_all) ** 2) * pdf_n_a_all)
                n_a_std_all = np.sqrt(n_a_var_all)
            else:
                pdf_n_a_all = hist_n_a_all
                n_a_mean_all = 0
                n_a_std_all = 0
            
            # n_pの統計情報を計算して保存
            if hist_n_p_all.sum() > 0:
                pdf_n_p_all = hist_n_p_all / hist_n_p_all.sum()
                n_p_mean_all = np.sum(bin_centers * pdf_n_p_all)
                n_p_var_all = np.sum(((bin_centers - n_p_mean_all) ** 2) * pdf_n_p_all)
                n_p_std_all = np.sqrt(n_p_var_all)
            else:
                pdf_n_p_all = hist_n_p_all
                n_p_mean_all = 0
                n_p_std_all = 0
            
            # ビン情報を両方のグループに保存
            n_a_group.create_dataset('bin_centers', data=bin_centers)
            n_p_group.create_dataset('bin_centers', data=bin_centers)
            
            # n_aの情報を保存
            n_a_group.create_dataset('pdf', data=pdf_n_a_all)
            n_a_group.create_dataset('hist', data=hist_n_a_all)
            n_a_group.attrs['mean'] = n_a_mean_all
            n_a_group.attrs['std'] = n_a_std_all
            n_a_group.attrs['avg_mean'] = avg_n_a_mean
            n_a_group.attrs['avg_std'] = avg_n_a_std
            n_a_group.attrs['total_count'] = total_count_n_a
            
            # n_pの情報を保存
            n_p_group.create_dataset('pdf', data=pdf_n_p_all)
            n_p_group.create_dataset('hist', data=hist_n_p_all)
            n_p_group.attrs['mean'] = n_p_mean_all
            n_p_group.attrs['std'] = n_p_std_all
            n_p_group.attrs['avg_mean'] = avg_n_p_mean
            n_p_group.attrs['avg_std'] = avg_n_p_std
            n_p_group.attrs['total_count'] = total_count_n_p
            
            # 共通のメタデータ
            threading_group.attrs['nchains'] = nchains
            threading_group.attrs['num_files'] = len(args.input)
            
            print(f"全体の統計情報を計算しました:")
            print(f"n_a 平均: {n_a_mean_all}, 標準偏差: {n_a_std_all}, 総数: {total_count_n_a}")
            print(f"n_p 平均: {n_p_mean_all}, 標準偏差: {n_p_std_all}, 総数: {total_count_n_p}")
            print(f"個別ファイルの平均の平均 n_a: {avg_n_a_mean}, 標準偏差の平均: {avg_n_a_std}")
            print(f"個別ファイルの平均の平均 n_p: {avg_n_p_mean}, 標準偏差の平均: {avg_n_p_std}")


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
