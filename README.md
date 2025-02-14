# HomologicalThreading
環状高分子のスレッディングをパーシステントホモロジーによって定量化するためのプログラム．

## Installation
### Requirements
- uv (Python プロジェクト管理ツール)
```bash
pip install uv
```
などでインストールしてください．

- Python 3.13 or later

- Fortran compiler (gfortran など)

ifx などの Fortran コンパイラが必要です.
デフォルトでは gfortran をMakefileで指定していますが，必要に応じて変更してください．

- CGAL 

システム全体にインストールする場合は単に，apt や brew などでインストールしてください．
```bash
sudo apt install libcgal-dev
```
管理者権限がない場合は，ユーザーローカルにインストールすることもできます．
いくつかの外部ライブラリが必要です．
    - BOOST

    基本的には公式HPのインストール方法に従ってください．
    https://www.boost.org/doc/libs/1_79_0/more/getting_started/unix-variants.html
    など．
    ```bash
    wget https://archives.boost.io/release/1.79.0/source/boost_1_79_0.tar.bz2
    tar xvf boost_1_79_0.tar.bz2
    cd boost_1_79_0
    ./bootstrap.sh 
    ./b2 headers
    ```
    などでインストールした後，`boost_1_79_0` があるディレクトリに $LD_LIBRARY_PATH を通してください．

    - CGAL
    こちらも公式HPのインストール方法に従ってください．
    ```bash
    wget https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2.tar.xz
    tar xvf CGAL-5.6.2.tar.xz
    ```
    などでインストールした後，`CGAL-5.6.2` があるディレクトリに $LD_LIBRARY_PATH を通しておくといいかも．

- Python 仮想環境

uv で仮想環境を作成してください．
```bash
uv sync
```
もし，HomCloud のビルドに失敗する場合は，多分CGALのライブラリが見つからないからだと思うので，
```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```
とするといいかもしれません．
