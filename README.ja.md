# Homological Threading ドキュメント

- [日本語版](README.ja.md)
- [English Version](README.md)

# Homological Threading

本プロジェクトは、**環状高分子**のスレッディング（絡み合い）を**パーシステントホモロジー**に基づいて定量化するためのプログラムです。パーシステントホモロジーにより、データの幾何学的・位相的な特徴（穴、空洞、連結成分など）を多段階で解析し、環状高分子間の複雑な絡み合いを捉えます。

---

## 目次
- [Homological Threading ドキュメント](#homological-threading-ドキュメント)
- [Homological Threading](#homological-threading)
  - [目次](#目次)
  - [特徴](#特徴)
  - [ディレクトリ構造](#ディレクトリ構造)
  - [3. 動作環境と依存ライブラリ](#3-動作環境と依存ライブラリ)
    - [3.1 必要条件](#31-必要条件)
      - [3.1.1 Python環境](#311-python環境)
      - [3.1.2 Fortranコンパイラ](#312-fortranコンパイラ)
      - [3.1.3 CGAL (Computational Geometry Algorithms Library)](#313-cgal-computational-geometry-algorithms-library)
  - [4. インストール方法](#4-インストール方法)
    - [4.1 必要なライブラリの準備](#41-必要なライブラリの準備)
    - [4.2 Python仮想環境のセットアップ](#42-python仮想環境のセットアップ)
    - [4.3 プロジェクトのビルド](#43-プロジェクトのビルド)
  - [5. 使用方法](#5-使用方法)
    - [5.1 パーシステント図の計算と解析](#51-パーシステント図の計算と解析)
    - [5.2 Betti数の計算](#52-betti数の計算)
    - [5.3 結果の可視化](#53-結果の可視化)
  - [6. テストの実行](#6-テストの実行)
  - [7. トラブルシューティング](#7-トラブルシューティング)
    - [7.1 インストール関連の問題](#71-インストール関連の問題)
    - [7.2 実行時の問題](#72-実行時の問題)
  - [8. 開発者向け情報](#8-開発者向け情報)
    - [8.1 コードの拡張](#81-コードの拡張)
    - [8.2 リポジトリ構成の変更](#82-リポジトリ構成の変更)

---

## 特徴
- **パーシステントホモロジー解析**: 点群データをフィルトレーションにより多段階で解析し、連結成分、穴、空洞などの位相的特徴を定量化。
- **スレッディング定量化**: 単一の環状高分子および複数環状高分子間の絡み合い（スレッディング）を、パーシステント図やBetti数計算を通じて評価。
- **高速計算**: 計算コアはFortranで実装されており、大規模データにも迅速に対応可能。

---

## ディレクトリ構造

プロジェクトの主要なファイル・ディレクトリは以下の通りです：

```
HomologicalThreading/          ← プロジェクトルート
├── .gitignore
├── .python-version
├── build.sh                   ← ビルドスクリプト
├── main.py                    ← ルートレベルのエントリーポイント（インターフェース）
├── pyproject.toml             ← プロジェクト設定ファイル
├── README.ja.md               ← この日本語README
├── README.md                  ← 英語版README
├── uv.lock                    ← 依存関係ロックファイル
├── data/                     ← サンプルデータ等
│   ├── N10M100.data          ← LAMMPSデータファイルのサンプル
│   └── (その他のデータ)
├── docs/                     ← ドキュメント関連
│   ├── build/
│   └── source/
│       ├── _static/
│       └── _templates/
├── scripts/                  ← 解析・可視化スクリプト
│   ├── __init__.py
│   ├── analysis.py           ← パーシステント図計算とBetti数解析
│   └── plot.py               ← 結果の可視化（プロット機能を統合）
├── src/                      ← ソースコード
│   └── homological_threading/
│       ├── __init__.py
│       ├── lammps_io.py      ← LAMMPSデータファイルの入出力
│       ├── main.py           ← コア機能の実装
│       └── fortran/          ← 高速計算のFortran実装
│           ├── __init__.py
│           ├── compute.f90   ← スレッディング、Betti数計算の高速部分
│           └── Makefile
└── tests/                    ← テストコード
    └── test.py
```

※ READMEのディレクトリ構造は、最新のファイル構成に基づいて更新されています。

---

## 3. 動作環境と依存ライブラリ

### 3.1 必要条件

#### 3.1.1 Python環境
- Python 3.13以上
- uv（Pythonプロジェクト管理ツール）

uvのインストール:
```bash
pip install uv
```

#### 3.1.2 Fortranコンパイラ
- gfortran（GNU Fortran Compiler）または ifx（Intel Fortran Compiler）

macOSの場合:
```bash
brew install gcc
```

#### 3.1.3 CGAL (Computational Geometry Algorithms Library)
CGALはHomCloudライブラリの依存ライブラリとして必要です。

macOSの場合:
```bash
brew install cgal
```

または、ユーザーローカルにインストールする場合は、BOOSTとCGALをソースから構築してください。

---

## 4. インストール方法

### 4.1 必要なライブラリの準備
上記の依存ライブラリ（Python, Fortranコンパイラ, CGALなど）をインストールしてください。

### 4.2 Python仮想環境のセットアップ

uvを使用して仮想環境を作成・同期します:

```bash
uv sync
```

※ CGALのインクルードパスが必要な場合は、以下のように指定します:
```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

### 4.3 プロジェクトのビルド

Fortranコードのビルドやその他必要なビルド処理は以下のコマンドで実行します:

```bash
./build.sh
```

※ Fortranコンパイラが見つからない場合は、`src/homological_threading/fortran/Makefile`内の`FC`変数を適切に設定してください。

---

## 5. 使用方法

### 5.1 パーシステント図の計算と解析
LAMMPSデータファイルから環状高分子のパーシステント図（PD）を計算し、スレッディングを解析します:

```bash
uv run main.py analysis pd -i data/N10M100.data -o output_directory
```

このコマンドは以下の処理を実行します:
1. LAMMPSデータファイルから環状高分子の座標を読み込み
2. 単一環状高分子のパーシステント図（PD_i）を計算
3. 複数環状高分子間のスレッディング解析（PD_i_cup_jの計算等）
4. 解析結果をHDF5形式で保存

### 5.2 Betti数の計算

保存されたHDF5ファイルからBetti数を計算します:

```bash
uv run main.py analysis betti -i output_directory/*.h5 -f output_directory/analysis.h5
```

### 5.3 結果の可視化

解析結果の可視化は以下のコマンドで実行します:

```bash
uv run main.py plot -i output_directory/analysis.h5
```

---

## 6. テストの実行

テストコードは以下のコマンドで実行できます:

```bash
uv run tests/test.py
```

---

## 7. トラブルシューティング

### 7.1 インストール関連の問題

- **HomCloudのビルド失敗**: CGALのインクルードパスが正しいか確認してください。
- **Fortranコンパイル失敗**: `src/homological_threading/fortran/Makefile`内の`FC`変数が正しいか確認。

### 7.2 実行時の問題

- **計算速度が遅い場合**:
  - 環境変数 `OMP_NUM_THREADS` を設定してスレッド数を増やす（例：`export OMP_NUM_THREADS=8`）。
  - Pythonのマルチプロセス処理の有効化を検討。

---

## 8. 開発者向け情報

### 8.1 コードの拡張
新機能を追加する場合は、以下のファイルを編集してください:
- コア解析ロジック: `src/homological_threading/main.py`
- 入出力処理: `src/homological_threading/lammps_io.py`
- 高速計算部分: `src/homological_threading/fortran/compute.f90` および `Makefile`
- テストコード: `tests/test.py`

### 8.2 リポジトリ構成の変更
ルートレベルのファイル（`main.py`, `build.sh`, `pyproject.toml` など）も必要に応じて更新してください。

