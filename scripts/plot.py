import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import pathlib


def get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # PD command
    pd_parser = subparsers.add_parser("pd", help="Plot persistence diagrams")
    pd_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files", dest="filename")
    pd_parser.add_argument("-o", "--output", default=None, help="Output image file")
    pd_parser.add_argument("--xlim", nargs=2, type=float, help="X-axis limits [min max]", default=[0, 50])
    pd_parser.add_argument("--ylim", nargs=2, type=float, help="Y-axis limits [min max]", default=[0, 50])
    
    # Betti command
    betti_parser = subparsers.add_parser("betti", help="Plot Betti numbers")
    betti_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files", dest="filename")
    betti_parser.add_argument("-o", "--output", default=None, help="Output image file")
    betti_parser.add_argument("--xlim", nargs=2, type=float, help="X-axis limits [min max]", default=[0, 10])
    betti_parser.add_argument("--ylim", nargs=2, type=float, help="Y-axis limits [min max]")
    
    # Num threading command
    num_threading_parser = subparsers.add_parser("num_threading", help="Plot number of threading")
    num_threading_parser.add_argument("-i", "--input", nargs="+", help="Input HDF5 files", dest="filename")
    num_threading_parser.add_argument("-o", "--output", default=None, help="Output image file")
    num_threading_parser.add_argument("--type", choices=["hist", "pdf"], default="pdf", 
                                     help="Plot histogram or probability density function")
    num_threading_parser.add_argument("--log", action="store_true", help="Use log scale for y-axis")
    num_threading_parser.add_argument("--xlim", nargs=2, type=float, help="X-axis limits [min max]")
    num_threading_parser.add_argument("--ylim", nargs=2, type=float, help="Y-axis limits [min max]")

    return parser.parse_args()


def betti(args):
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 1, figsize=(10, 5*num_files))
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot_betti(axes[i], filename)
        axes[i].set_title(pathlib.Path(filename).stem)
        
        # Set axis limits if specified
        if args.xlim:
            axes[i].set_xlim(args.xlim)
        if args.ylim:
            axes[i].set_ylim(args.ylim)
    
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output, dpi=300)
    else:
        plt.show()


def plot_betti(ax, filename):
    # HDF5からベッチ数を読み込む
    with h5py.File(filename, "r") as f:
        if 'betti' in f:
            betti_group = f['betti']
            alphas = betti_group['alphas'][:]
            betti_pd_i = betti_group['pd_i'][:]
            betti_threading = betti_group['threading'][:]
            
            # オプションでpd_i_cup_jも表示
            # if 'pd_i_cup_j' in betti_group:
            #     betti_pd_i_cup_j = betti_group['pd_i_cup_j'][:]
        else:
            # 古い形式のnpzファイルからの読み込みサポート
            data = np.load(filename)
            alphas = data["alphas"]
            betti_pd_i = data["betti_pd_i"]
            betti_threading = data["betti_threading"]
            if "betti_pd_i_cup_j" in data:
                betti_pd_i_cup_j = data["betti_pd_i_cup_j"]
    
    # プロット
    const_pd_i = 100.0  # ベッチ数のスケール調整定数
    ax.plot(np.sqrt(alphas), betti_pd_i / const_pd_i, label="i", color='blue')
    
    const_threading = 1.0  # threadingのスケール調整定数
    ax.plot(np.sqrt(alphas), betti_threading / const_threading, label="threading", color='red')
    
    # オプション：pd_i_cup_jも表示
    try:
        const_pd_i_cup_j = 100.0
        ax.plot(np.sqrt(alphas), betti_pd_i_cup_j / const_pd_i_cup_j, label="i cup j", color='green', alpha=0.7)
    except:
        pass

    ax.set_xlim(0, 10)
    ax.set_xlabel(r"$\sqrt{\alpha}$")
    ax.set_ylabel(r"$\beta(\alpha)$")
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()


def num_threading(args):
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 2, figsize=(15, 5*num_files))
    if num_files == 1:
        axes = np.array([axes])
    
    for i, filename in enumerate(args.filename):
        plot_num_threading(axes[i, 0], axes[i, 1], filename, plot_type=args.type, use_log=args.log, xlim=args.xlim, ylim=args.ylim)
        axes[i, 0].set_title(f"{pathlib.Path(filename).stem} - n_a")
        axes[i, 1].set_title(f"{pathlib.Path(filename).stem} - n_p")
    
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output, dpi=300)
    else:
        plt.show()


def plot_num_threading(ax_n_a, ax_n_p, filename, plot_type="pdf", use_log=False, xlim=None, ylim=None):
    with h5py.File(filename, "r") as f:
        if 'num_threading' not in f:
            print(f"Error: No num_threading data found in {filename}")
            return
        
        threading_group = f['num_threading']
        
        # n_aのデータを取得
        if 'n_a' in threading_group:
            n_a_group = threading_group['n_a']
            bin_centers_n_a = n_a_group['bin_centers'][:]
            
            if plot_type == "pdf":
                data_n_a = n_a_group['pdf'][:]
                ylabel = "Probability Density"
                # Use line plot for PDF
                ax_n_a.plot(bin_centers_n_a, data_n_a, '-', color='blue', linewidth=2, marker='o', markersize=4, alpha=0.7)
            else:  # hist
                data_n_a = n_a_group['hist'][:]
                ylabel = "Frequency"
                # Use bar plot for histogram
                ax_n_a.bar(bin_centers_n_a, data_n_a, width=0.8, alpha=0.7, color='blue')
            
            # 統計情報
            mean_n_a = n_a_group.attrs['mean']
            std_n_a = n_a_group.attrs['std']
            
            # Add vertical line for mean
            ax_n_a.axvline(mean_n_a, color='red', linestyle='--', 
                         label=f'Mean: {mean_n_a:.2f}, Std Dev: {std_n_a:.2f}')
            if use_log:
                ax_n_a.set_yscale('log')
            ax_n_a.set_xlabel('n_a (Number of Active threadings)')
            ax_n_a.set_ylabel(ylabel)
            ax_n_a.grid(True, linestyle='--', alpha=0.5)
            ax_n_a.legend()
            
            # Set axis limits if specified
            if xlim:
                ax_n_a.set_xlim(xlim)
            if ylim:
                ax_n_a.set_ylim(ylim)
        
        # n_pのデータを取得
        if 'n_p' in threading_group:
            n_p_group = threading_group['n_p']
            bin_centers_n_p = n_p_group['bin_centers'][:]
            
            if plot_type == "pdf":
                data_n_p = n_p_group['pdf'][:]
                ylabel = "Probability Density"
                # Use line plot for PDF
                ax_n_p.plot(bin_centers_n_p, data_n_p, '-', color='green', linewidth=2, marker='o', markersize=4, alpha=0.7)
            else:  # hist
                data_n_p = n_p_group['hist'][:]
                ylabel = "Frequency"
                # Use bar plot for histogram
                ax_n_p.bar(bin_centers_n_p, data_n_p, width=0.8, alpha=0.7, color='green')
            
            # 統計情報
            mean_n_p = n_p_group.attrs['mean']
            std_n_p = n_p_group.attrs['std']
            
            # Add vertical line for mean
            ax_n_p.axvline(mean_n_p, color='red', linestyle='--', 
                         label=f'Mean: {mean_n_p:.2f}, Std Dev: {std_n_p:.2f}')
            if use_log:
                ax_n_p.set_yscale('log')
            ax_n_p.set_xlabel('n_p (Number of Passive threadings)')
            ax_n_p.set_ylabel(ylabel)
            ax_n_p.grid(True, linestyle='--', alpha=0.5)
            ax_n_p.legend()
            
            # Set axis limits if specified
            if xlim:
                ax_n_p.set_xlim(xlim)
            if ylim:
                ax_n_p.set_ylim(ylim)


def pd(args):
    num_files = len(args.filename)
    fig, axes = plt.subplots(num_files, 3, figsize=(15, 5*num_files))
    if num_files == 1:
        axes = [axes]
    for i, filename in enumerate(args.filename):
        plot_pd(axes[i], filename)
    
    # Set axis limits
    for ax in axes:
        for a in ax:
            a.set_xlim(args.xlim)
            a.set_ylim(args.ylim)
    
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output, dpi=300)
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
        
        # 対角線を追加
        for a in ax:
            lims = [
                np.min([a.get_xlim(), a.get_ylim()]),
                np.max([a.get_xlim(), a.get_ylim()]),
            ]
            a.plot(lims, lims, 'k-', alpha=0.3, zorder=0)
            a.set_xlabel('Birth')
            a.set_ylabel('Death')


def main(args):
    if args.command == "pd":
        pd(args)
    elif args.command == "betti":
        betti(args)
    elif args.command == "num_threading":
        num_threading(args)
    else:
        raise ValueError("Unknown command")


if __name__ == "__main__":
    args = get_args()
    main(args)
