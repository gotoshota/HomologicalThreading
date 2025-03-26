# HomologicalThreading

環状高分子のスレッディングをパーシステントホモロジーによって定量化するためのプログラム．

## 1. プロジェクト概要

### 1.1 背景

環状高分子（リングポリマー）は、一般的な線状高分子と異なり、分子鎖の両端が結合して環状構造を形成しています。これらの環状高分子が互いに「スレッディング」する現象（一方の高分子が他方の高分子の輪を通り抜ける状態）は、物質の物理的性質に大きな影響を与えます。

しかし、このスレッディング現象は従来の方法では定量的に検出することが困難でした。本プロジェクトでは、トポロジカルデータ解析の一手法である「パーシステントホモロジー」を用いて、このスレッディング現象を数学的に検出・定量化する手法を実装しています。

### 1.2 パーシステントホモロジーとは

パーシステントホモロジーは、データの位相的特徴（穴や空洞など）を異なるスケールで捉える数学的手法です。点群データに対して、点間の距離を徐々に大きくしていき（フィルトレーション）、その過程で現れる・消える位相的特徴を追跡します。

この手法により、環状高分子の形状とその相互作用を数学的に記述し、スレッディング現象を定量的に評価することが可能になります。

### 1.3 本プロジェクトの目的

本プロジェクトの主な目的は以下の通りです：

1. 環状高分子の構造をパーシステントホモロジーで解析する
2. 環状高分子間のスレッディングを数学的に検出する
3. スレッディングの度合いを定量化する

## 2. プロジェクト構造

### 2.1 ディレクトリ構造

```
HomologicalThreading/
├── data/                  # サンプルデータ
│   ├── ht.h5              # 解析結果のサンプル
│   └── N10M100.data       # LAMMPSデータファイルのサンプル
├── scripts/               # 解析・可視化スクリプト
│   ├── analysis.py        # パーシステント図計算とベッティ数解析
│   ├── plot_betti.py      # ベッティ数のプロット
│   └── plot_pd.py         # パーシステント図の可視化
├── src/                   # ソースコード
│   └── homological_threading/
│       ├── __init__.py
│       ├── lammps_io.py   # LAMMPSデータファイルの入出力
│       ├── main.py        # メインの実装
│       └── fortran/       # Fortranによる高速化実装
│           ├── __init__.py
│           ├── compute.f90 # 計算コア部分
│           └── Makefile
├── tests/                 # テストコード
│   └── test.py
├── .gitignore
├── .python-version
├── build.sh               # ビルドスクリプト
├── pyproject.toml         # Pythonプロジェクト設定
├── README.md              # このファイル
└── uv.lock                # 依存関係ロックファイル
```

### 2.2 主要なコンポーネント

#### 2.2.1 Python部分

- `homological_threading/main.py`: 
  - `HomologicalThreading`クラス: プロジェクトの中心となるクラス
  - `PD_i`クラス: 単一環状高分子のパーシステント図を計算
  - `PD_i_cup_j`クラス: 2つの環状高分子のカップ積のパーシステント図を計算
  - `Threading`クラス: スレッディングの検出と定量化

- `homological_threading/lammps_io.py`: 
  - `LammpsData`クラス: LAMMPSデータファイルの読み書き
  - `polyWrap`メソッド: 周期境界条件での分子の適切な配置

#### 2.2.2 Fortran部分

- `homological_threading/fortran/compute.f90`: 
  - `threading`サブルーチン: スレッディング計算の高速実装
  - `betti_number`サブルーチン: ベッティ数計算の高速実装

#### 2.2.3 スクリプト

- `scripts/analysis.py`: パーシステント図の計算とベッティ数の解析
- `scripts/plot_pd.py`: パーシステント図の可視化
- `scripts/plot_betti.py`: ベッティ数のプロット

## 3. インストール方法

### 3.1 必要条件

#### 3.1.1 Python環境

- Python 3.13以上
- uv (Pythonプロジェクト管理ツール)

uvのインストール:
```bash
pip install uv
```

#### 3.1.2 Fortranコンパイラ

以下のいずれかのFortranコンパイラが必要です:
- gfortran (GNU Fortran Compiler)
- ifx (Intel Fortran Compiler)

Ubuntuの場合:
```bash
sudo apt install gfortran
```

macOSの場合:
```bash
brew install gcc
```

#### 3.1.3 CGAL (Computational Geometry Algorithms Library)

HomCloudライブラリの依存ライブラリとして必要です。

**システム全体にインストールする場合:**

Ubuntuの場合:
```bash
sudo apt install libcgal-dev
```

macOSの場合:
```bash
brew install cgal
```

**ユーザーローカルにインストールする場合:**

管理者権限がない場合は、以下の手順でユーザーローカルにインストールできます。

1. BOOSTのインストール:
   ```bash
   wget https://archives.boost.io/release/1.79.0/source/boost_1_79_0.tar.bz2
   tar xvf boost_1_79_0.tar.bz2
   cd boost_1_79_0
   ./bootstrap.sh 
   ./b2 headers
   ```
   インストール後、環境変数を設定:
   ```bash
   export LD_LIBRARY_PATH=/path/to/boost_1_79_0:$LD_LIBRARY_PATH
   ```

2. CGALのインストール:
   ```bash
   wget https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2.tar.xz
   tar xvf CGAL-5.6.2.tar.xz
   ```
   インストール後、環境変数を設定:
   ```bash
   export LD_LIBRARY_PATH=/path/to/CGAL-5.6.2:$LD_LIBRARY_PATH
   ```

### 3.2 Python仮想環境のセットアップ

uvを使用して仮想環境を作成します:

```bash
uv sync
```

HomCloudのビルドに失敗する場合は、CGALのインクルードパスを指定してみてください:

```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

### 3.3 ビルド

以下のコマンドでプロジェクトをビルドします:

```bash
./build.sh
```

Fortranコンパイラが見つからない場合は、`src/homological_threading/fortran/Makefile`の`FC`変数を適切なコンパイラに設定してください。

例:
```makefile
# gfortranを使用する場合
FC = gfortran

# Intel Fortranコンパイラを使用する場合
# FC = ifx
```

## 4. 使用方法

### 4.1 基本的な使用例

#### 4.1.1 パーシステント図の計算

LAMMPSデータファイルからパーシステント図を計算します:

```bash
python scripts/analysis.py pd -i data/N10M100.data -o output_directory
```

このコマンドは以下の処理を行います:
1. LAMMPSデータファイルから環状高分子の座標を読み込む
2. 単一環状高分子のパーシステント図（PD_i）を計算
3. 環状高分子ペアのパーシステント図（PD_i_cup_j）を計算
4. スレッディングの検出と定量化
5. 結果をHDF5ファイルに保存

#### 4.1.2 ベッティ数の計算

保存されたHDF5ファイルからベッティ数を計算します:

```bash
python scripts/analysis.py betti -i output_directory/*.h5 -o output_directory
```

#### 4.1.3 結果の可視化

パーシステント図の可視化:

```bash
python scripts/plot_pd.py output_directory/*.h5
```

ベッティ数のプロット:

```bash
python scripts/plot_betti.py output_directory/betti.h5
```

### 4.2 入力データ形式

本プロジェクトはLAMMPSデータファイル形式の入力を受け付けます。このファイルには、環状高分子の原子座標、結合情報、周期境界条件などが含まれています。

サンプルデータとして`data/N10M100.data`が提供されています。このファイルは10個のビーズからなる環状高分子が100個含まれるシステムを表しています。

### 4.3 出力データの解釈

#### 4.3.1 HDF5ファイル構造

出力されるHDF5ファイルには以下の情報が含まれています:

- `/pd_i/pd`: 単一環状高分子のパーシステント図
- `/pd_i_cup_j/pd`: 環状高分子ペアのパーシステント図
- `/threading/flags`: スレッディングの有無を示すフラグ
- `/threading/pd`: スレッディングに関連するパーシステント図
- `/Metadata`: 解析に関するメタデータ

#### 4.3.2 パーシステント図の解釈

パーシステント図は、位相的特徴の「誕生」と「消滅」のスケールを表します。横軸が誕生スケール、縦軸が消滅スケールです。対角線から離れた点ほど、「持続性の高い」特徴を表します。

#### 4.3.3 スレッディングフラグの解釈

`threading/flags`は、環状高分子ペア間のスレッディングの有無を示す行列です。`flags[i, j] = True`は、高分子jが高分子iをスレッディングしていることを示します。

## 5. 理論的背景

### 5.1 パーシステントホモロジーの基礎

パーシステントホモロジーは、データの位相的特徴を異なるスケールで捉える手法です。点群に対して、点間の距離εを徐々に大きくしていき、その過程で形成される単体複体（シンプリシャル・コンプレックス）の位相的特徴を追跡します。

- 0次元ホモロジー: 連結成分（点）
- 1次元ホモロジー: 輪（穴）
- 2次元ホモロジー: 空洞

環状高分子の場合、1次元ホモロジーが特に重要です。

### 5.2 スレッディングの検出方法

本プロジェクトでは、以下の手順でスレッディングを検出します:

1. 単一環状高分子iのパーシステント図（PD_i）を計算
2. 環状高分子iとjのカップ積のパーシステント図（PD_i_cup_j）を計算
3. PD_iとPD_i_cup_jの差分を計算
4. 差分が存在する場合、高分子jが高分子iをスレッディングしていると判定

この方法は、スレッディングによって生じる位相的変化を捉えることができます。

## 6. トラブルシューティング

### 6.1 インストール関連の問題

#### 6.1.1 HomCloudのビルドに失敗する場合

CGALのライブラリが見つからない可能性があります。以下を試してください:

```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

#### 6.1.2 Fortranコンパイルに失敗する場合

`src/homological_threading/fortran/Makefile`の`FC`変数が適切に設定されているか確認してください。

### 6.2 実行時の問題

#### 6.2.1 メモリ不足エラー

大規模なシステムを解析する場合、メモリ不足になる可能性があります。以下の対策を試してください:

- 小さなサブセットで解析を行う
- マルチプロセス処理を無効にする（`mp=False`オプションを使用）

#### 6.2.2 計算速度が遅い場合

- OpenMPのスレッド数を増やす: `export OMP_NUM_THREADS=8`
- Pythonのマルチプロセス処理を有効にする（`mp=True`オプションを使用）

## 7. 開発者向け情報

### 7.1 コードの拡張

新しい機能を追加する場合は、以下のファイルを編集してください:

- 新しい解析手法: `src/homological_threading/main.py`
- 新しい入出力形式: `src/homological_threading/lammps_io.py`
- 計算の高速化: `src/homological_threading/fortran/compute.f90`

### 7.2 テスト

テストを実行するには:

```bash
python tests/test.py
```

## 8. 参考文献

1. パーシステントホモロジーの理論:
   - Edelsbrunner, H., & Harer, J. (2010). Computational Topology: An Introduction. American Mathematical Society.

2. 先行研究:
   - F. Landuzzi, T. Nakamura, D. Michieletto, and T. Sakaue, Persistent homology of entangled rings, Physical Review Research 2, 033529 (2020).

3. HomCloudライブラリ:
   - Obayashi, I. (2018). HomCloud: A computing software for persistent homology. https://homcloud.dev/
