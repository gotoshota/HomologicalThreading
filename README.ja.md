# Homological Threading Documentation

- [日本語版](README.ja.md)
- [English Version](README.md)

# Homological Threading

このプロジェクトは、**環状高分子**のスレッディングを**パーシステントホモロジー**に基づいて定量化するためのプログラムです。  
パーシステントホモロジーを用いることで、データの幾何学的および位相的特徴（穴、空洞、連結成分など）を異なるスケールで解析でき、環状高分子間の複雑な絡み合い（スレッディング）を捉えることが可能です。

---

## 目次
- [特徴](#特徴)
- [ディレクトリ構造](#ディレクトリ構造)
- [動作環境と依存ライブラリ](#動作環境と依存ライブラリ)
- [インストール方法](#インストール方法)
  - [1. 必要なライブラリの準備](#1-必要な-ライブラリの-準備)
  - [2. Python仮想環境のセットアップ](#2-python仮想環境の-セットアップ)
  - [3. プロジェクトのビルド](#3-プロジェクトの-ビルド)
- [使用方法](#使用方法)
  - [パーシステント図の計算と解析](#パーシステント図の-計算と解析)
  - [Betti数の計算](#Betti数の-計算)
  - [結果の可視化](#結果の-可視化)
- [テストの実行](#テストの-実行)
- [トラブルシューティング](#トラブルシューティング)
- [貢献について](#貢献について)
- [ライセンス](#ライセンス)

---

## 特徴
- **パーシステントホモロジー解析**: フィルトレーションを通じて点群データの位相的特徴（穴、連結成分、空洞など）を追跡し、定量的な解析を可能にします。
- **スレッディング定量化**: 単一の環状高分子および複数の環状高分子間の絡み合いを、パーシステント図（PD）やbetti数計算を通じて評価します。
- **高速計算**: 計算コアはFortranにより実装され、大規模データに対しても効率的に演算可能です。

---

## ディレクトリ構造
以下はプロジェクトの主なディレクトリ構造です：

```
HomologicalThreading/
├── data/                  # サンプルデータ
│   ├── ht.h5              # 解析結果のサンプル
│   └── N10M100.data       # LAMMPSデータファイルのサンプル
├── scripts/               # 解析・可視化スクリプト
│   ├── analysis.py        # パーシステント図計算とbetti数解析
│   ├── plot_betti.py      # betti数のプロット
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
├── main.py                # インターフェース
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
  - `betti_number`サブルーチン: betti数計算の高速実装

#### 2.2.3 スクリプト

- `main.py`: プログラム全体のインターフェースとして機能し、ユーザの入力に応じて解析やプロットの各サブモジュールを呼び出すエントリーポイント
- `scripts/analysis.py`: パーシステント図の計算とbetti数の解析
- `scripts/plot_pd.py`: パーシステント図の可視化
- `scripts/plot_betti.py`: betti数のプロット

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

Fortranコンパイラが必要です:
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
uv run python main.py analysis pd -i data/N10M100.data -o output_directory
```

このコマンドは以下の処理を行います:
1. LAMMPSデータファイルから環状高分子の座標を読み込む
2. 単一環状高分子のパーシステント図（PD_i）を計算
3. 環状高分子ペアのパーシステント図（PD_i_cup_j）を計算
4. スレッディングの検出と定量化
5. 結果をHDF5ファイルに保存

#### 4.1.2 betti数の計算

保存されたHDF5ファイルからbetti数を計算します:

```bash
uv run main analysis betti -i output_directory/*.h5 -f output_directory/analysis.h5
```

#### 4.1.3 結果の可視化

パーシステント図の可視化:

```bash
uv run main.py plot pd -i output_directory/analysis.h5
uv run main.py plot pd -i output_directory/*.h5
```

betti数のプロット:

```bash
uv run main plot betti -i data/analysis.h5
uv run main plot betti -i data/*.h5
```

### 4.2 入力データ形式

本プロジェクトはLAMMPSデータファイル形式の入力を受け付けます。このファイルには、環状高分子の原子座標、結合情報、周期境界条件などが含まれています。

サンプルデータとして`data/N10M100.data`が提供されています。このファイルは10個のビーズからなる環状高分子が100個含まれるシステムを表しています。

### 4.3 出力データの解釈

#### 4.3.1 HDF5ファイル構造

出力されるHDF5ファイルには以下の情報が含まれています:

- `/pd_i/pd`: 単一環状高分子のパーシステント図
    - shape: [num_chains, num_points, 2]
        - num_chains: 環状高分子の数
        - num_points: 点の数
        - 2: birth, deathの2次元座標
    - 例1: `pd_i/pd[0]`は1番目の環状高分子のパーシステント図
    - 例2: `pd_i/pd[0, 0]`は1番目の環状高分子の1番目の点のパーシステント図
    - 例3: `pd_i/pd[0, 0, 0]`は1番目の環状高分子の1番目の点のbirth, `pd_i/pd[0, 0, 1]`は1番目の環状高分子の1番目の点のdeath
- `/pd_i_cup_j/pd`: 環状高分子ペアのパーシステント図
    - shape: [num_chains, num_chains, num_points, 2]
        - num_chains: 環状高分子の数
        - num_points: 点の数
        - 2: birth, deathの2次元座標
    - 例1: `pd_i_cup_j/pd[0, 1]`は1番目の環状高分子と2番目の環状高分子のパーシステント図
    - 例2: `pd_i_cup_j/pd[0, 1, 0]`は1番目の環状高分子と2番目の環状高分子の1番目の点のパーシステント図
    - 例3: `pd_i_cup_j/pd[0, 1, 0, 0]`は1番目の環状高分子と2番目の環状高分子の1番目の点のbirth, `pd_i_cup_j/pd[0, 1, 0, 1]`は1番目の環状高分子と2番目の環状高分子の1番目の点のdeath
- `/threading/flags`: スレッディングの有無を示すフラグ
    - shape: [num_chains, num_chains]
        - num_chains: 環状高分子の数
            - True, False で構成されるboolean行列
            - 順に， passive, active をあらわす
    - 例: `threading/flags[0, 1]`は2番目の環状高分子が1番目の環状高分子をスレッディングしているかどうかを示す
- `/threading/pd`: スレッディングに関連するパーシステント図
    - shape: [num_chains, num_chains, num_points, 2]
        - num_chains: 環状高分子の数
        - num_points: 点の数
        - 2: birth, deathの2次元座標
    - 例1: `threading/pd[0, 1]`は1番目の環状高分子と2番目の環状高分子のスレッディングに関連するパーシステント図
- `/Metadata`: 解析に関するメタデータ

#### 4.3.2 パーシステント図の解釈

パーシステント図は、位相的特徴の「誕生」と「消滅」のスケールを表します。横軸が誕生スケール、縦軸が消滅スケールです。対角線から離れた点ほど、半径パラメータの変化に対して「持続性の高い」特徴を表します。

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

CGALのインクルードパスが正しく設定されていない場合、HomCloudのビルドに失敗することがあります。
正しいパスが設定されているか確認してください。
例えば，以下を試してください:

```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

#### 6.1.2 Fortranコンパイルに失敗する場合

`src/homological_threading/fortran/Makefile`の`FC`変数が適切に設定されているか確認してください。

### 6.2 実行時の問題

#### 6.2.1 計算速度が遅い場合

- OpenMPのスレッド数を増やす: `export OMP_NUM_THREADS=8`
- Pythonのマルチプロセス処理を有効にする（`mp=True`オプションを使用）

## 7. 開発者向け情報

### 7.1 コードの拡張

新しい機能を追加する場合は、以下のファイルを編集してください:

- 新しい解析: `src/homological_threading/main.py`
- 入出力形式の修正: `src/homological_threading/lammps_io.py`
- 計算の高速化: `src/homological_threading/fortran/compute.f90`

### 7.2 テスト

テストを実行するには:

```bash
python tests/test.py
