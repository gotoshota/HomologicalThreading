# Homological Threading Documentation

- [日本語版](README.ja.md)
- [English Version](README.md)

# Homological Threading

This project is a program to quantify the threading of **cyclic polymers** based on **persistent homology**.  
By using persistent homology, it is possible to analyze the geometric and topological features (holes, cavities, connected components, etc.) of data at different scales, capturing the complex entanglement (threading) between cyclic polymers.

---

## Table of Contents
- [Features](#features)
- [Directory Structure](#directory-structure)
- [Environment and Dependencies](#environment-and-dependencies)
- [Installation](#installation)
  - [1. Prepare Required Libraries](#1-prepare-required-libraries)
  - [2. Set Up Python Virtual Environment](#2-set-up-python-virtual-environment)
  - [3. Build the Project](#3-build-the-project)
- [Usage](#usage)
  - [Calculate and Analyze Persistent Diagrams](#calculate-and-analyze-persistent-diagrams)
  - [Calculate Betti Numbers](#calculate-betti-numbers)
  - [Visualize Results](#visualize-results)
- [Run Tests](#run-tests)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

---

## Features
- **Persistent Homology Analysis**: Tracks topological features (holes, connected components, cavities, etc.) of point cloud data through filtration, enabling quantitative analysis.
- **Threading Quantification**: Evaluates the entanglement between single and multiple cyclic polymers through persistent diagrams (PD) and Betti number calculations.
- **High-Speed Computation**: The computation core is implemented in Fortran, allowing efficient operations even on large-scale data.

---

## Directory Structure
Below is the main directory structure of the project:

```
HomologicalThreading/
├── data/                  # Sample data
│   ├── ht.h5              # Sample analysis result
│   └── N10M100.data       # Sample LAMMPS data file
├── scripts/               # Analysis and visualization scripts
│   ├── analysis.py        # Persistent diagram calculation and Betti number analysis
│   ├── plot_betti.py      # Betti number plotting
│   └── plot_pd.py         # Persistent diagram visualization
├── src/                   # Source code
│   └── homological_threading/
│       ├── __init__.py
│       ├── lammps_io.py   # LAMMPS data file I/O
│       ├── main.py        # Main implementation
│       └── fortran/       # Fortran-based high-speed implementation
│           ├── __init__.py
│           ├── compute.f90 # Computation core
│           └── Makefile
├── tests/                 # Test code
│   └── test.py
├── .gitignore
├── .python-version
├── build.sh               # Build script
├── pyproject.toml         # Python project settings
├── README.md              # This file
├── main.py                # Interface
└── uv.lock                # Dependency lock file

```

### 2.2 Main Components

#### 2.2.1 Python Part

- `homological_threading/main.py`: 
  - `HomologicalThreading` class: Central class of the project
  - `PD_i` class: Calculates the persistent diagram of a single cyclic polymer
  - `PD_i_cup_j` class: Calculates the persistent diagram of the cup product of two cyclic polymers
  - `Threading` class: Detection and quantification of threading

- `homological_threading/lammps_io.py`: 
  - `LammpsData` class: Reading and writing LAMMPS data files
  - `polyWrap` method: Proper placement of molecules under periodic boundary conditions

#### 2.2.2 Fortran Part

- `homological_threading/fortran/compute.f90`: 
  - `threading` subroutine: High-speed implementation of threading calculation
  - `betti_number` subroutine: High-speed implementation of Betti number calculation

#### 2.2.3 Scripts

- `main.py`: Functions as the interface for the entire program, calling each submodule for analysis and plotting based on user input
- `scripts/analysis.py`: Calculation of persistent diagrams and Betti number analysis
- `scripts/plot_pd.py`: Visualization of persistent diagrams
- `scripts/plot_betti.py`: Plotting of Betti numbers

## 3. Installation

### 3.1 Requirements

#### 3.1.1 Python Environment

- Python 3.13 or higher
- uv (Python project management tool)

Install uv:
```bash
pip install uv
```

#### 3.1.2 Fortran Compiler

A Fortran compiler is required:
- gfortran (GNU Fortran Compiler)
- ifx (Intel Fortran Compiler)

On Ubuntu:
```bash
sudo apt install gfortran
```

On macOS:
```bash
brew install gcc
```

#### 3.1.3 CGAL (Computational Geometry Algorithms Library)

Required as a dependency for the HomCloud library.

**To install system-wide:**

On Ubuntu:
```bash
sudo apt install libcgal-dev
```

On macOS:
```bash
brew install cgal
```

**To install locally:**

If you do not have administrator privileges, you can install locally using the following steps.

1. Install BOOST:
   ```bash
   wget https://archives.boost.io/release/1.79.0/source/boost_1_79_0.tar.bz2
   tar xvf boost_1_79_0.tar.bz2
   cd boost_1_79_0
   ./bootstrap.sh 
   ./b2 headers
   ```
   After installation, set the environment variable:
   ```bash
   export LD_LIBRARY_PATH=/path/to/boost_1_79_0:$LD_LIBRARY_PATH
   ```

2. Install CGAL:
   ```bash
   wget https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2.tar.xz
   tar xvf CGAL-5.6.2.tar.xz
   ```
   After installation, set the environment variable:
   ```bash
   export LD_LIBRARY_PATH=/path/to/CGAL-5.6.2:$LD_LIBRARY_PATH
   ```

### 3.2 Set Up Python Virtual Environment

Create a virtual environment using uv:

```bash
uv sync
```

If building HomCloud fails, try specifying the CGAL include path:

```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

### 3.3 Build

Build the project with the following command:

```bash
./build.sh
```

If the Fortran compiler is not found, set the `FC` variable in `src/homological_threading/fortran/Makefile` to the appropriate compiler.

Example:
```makefile
# Use gfortran
FC = gfortran

# Use Intel Fortran Compiler
# FC = ifx
```

## 4. Usage

### 4.1 Basic Usage Examples

#### 4.1.1 Calculate Persistent Diagrams

Calculate persistent diagrams from LAMMPS data files:

```bash
python scripts/analysis.py pd -i data/N10M100.data -o output_directory
```

This command performs the following:
1. Reads coordinates of cyclic polymers from LAMMPS data files
2. Calculates the persistent diagram (PD_i) of a single cyclic polymer
3. Calculates the persistent diagram (PD_i_cup_j) of cyclic polymer pairs
4. Detects and quantifies threading
5. Saves results to an HDF5 file

#### 4.1.2 Calculate Betti Numbers

Calculate Betti numbers from saved HDF5 files:

```bash
python main analysis betti -i output_directory/*.h5 -f output_directory/analysis.h5
```

#### 4.1.3 Visualize Results

Visualize persistent diagrams:

```bash
uv run scripts/plot.py pd -i output_directory/analysis.h5
uv run scripts/plot.py pd -i output_directory/*.h5
```

Plot Betti numbers:

```bash
uv run main plot betti -i data/analysis.h5
uv run main plot betti -i data/*.h5
```

### 4.2 Input Data Format

This project accepts input in the form of LAMMPS data files. These files contain atomic coordinates, bonding information, and periodic boundary conditions for cyclic polymers.

Sample data `data/N10M100.data` is provided, representing a system with 100 cyclic polymers each consisting of 10 beads.

### 4.3 Interpretation of Output Data

#### 4.3.1 HDF5 File Structure

The output HDF5 file contains the following information:

- `/pd_i/pd`: Persistent diagram of a single cyclic polymer
    - shape: [num_chains, num_points, 2]
        - num_chains: Number of cyclic polymers
        - num_points: Number of points
        - 2: 2D coordinates of birth and death
    - Example 1: `pd_i/pd[0]` is the persistent diagram of the first cyclic polymer
    - Example 2: `pd_i/pd[0, 0]` is the persistent diagram of the first point of the first cyclic polymer
    - Example 3: `pd_i/pd[0, 0, 0]` is the birth of the first point of the first cyclic polymer, `pd_i/pd[0, 0, 1]` is the death of the first point of the first cyclic polymer
- `/pd_i_cup_j/pd`: Persistent diagram of cyclic polymer pairs
    - shape: [num_chains, num_chains, num_points, 2]
        - num_chains: Number of cyclic polymers
        - num_points: Number of points
        - 2: 2D coordinates of birth and death
    - Example 1: `pd_i_cup_j/pd[0, 1]` is the persistent diagram of the first and second cyclic polymers
    - Example 2: `pd_i_cup_j/pd[0, 1, 0]` is the persistent diagram of the first point of the first and second cyclic polymers
    - Example 3: `pd_i_cup_j/pd[0, 1, 0, 0]` is the birth of the first point of the first and second cyclic polymers, `pd_i_cup_j/pd[0, 1, 0, 1]` is the death of the first point of the first and second cyclic polymers
- `/threading/flags`: Flags indicating the presence of threading
    - shape: [num_chains, num_chains]
        - num_chains: Number of cyclic polymers
            - Boolean matrix composed of True and False
            - Represents passive and active, respectively
    - Example: `threading/flags[0, 1]` indicates whether the second cyclic polymer is threading the first cyclic polymer
- `/threading/pd`: Persistent diagram related to threading
    - shape: [num_chains, num_chains, num_points, 2]
        - num_chains: Number of cyclic polymers
        - num_points: Number of points
        - 2: 2D coordinates of birth and death
    - Example 1: `threading/pd[0, 1]` is the persistent diagram related to threading of the first and second cyclic polymers
- `/Metadata`: Metadata related to the analysis

#### 4.3.2 Interpretation of Persistent Diagrams

Persistent diagrams represent the "birth" and "death" scales of topological features. The horizontal axis represents the birth scale, and the vertical axis represents the death scale. Points farther from the diagonal represent features with high "persistence" against changes in the radius parameter.

#### 4.4 Reading HDF5 Files with Command Line Tool

Typically,

### 5.1 Basics of Persistent Homology

Persistent homology is a method for capturing the topological features of data at different scales. For point clouds, the distance ε between points is gradually increased, and the topological features of the simplicial complex formed during this process are tracked.

- 0-dimensional homology: Connected components (points)
- 1-dimensional homology: Loops (holes)
- 2-dimensional homology: Cavities

For cyclic polymers, 1-dimensional homology is particularly important.

### 5.2 Method for Detecting Threading

In this project, threading is detected using the following steps:

1. Calculate the persistent diagram (PD_i) of a single cyclic polymer i
2. Calculate the persistent diagram (PD_i_cup_j) of the cup product of cyclic polymers i and j
3. Calculate the difference between PD_i and PD_i_cup_j
4. If a difference exists, it is determined that polymer j is threading polymer i

This method captures the topological changes caused by threading.

## 6. Troubleshooting

### 6.1 Installation Issues

#### 6.1.1 If Building HomCloud Fails

If the CGAL include path is not set correctly, building HomCloud may fail.
Please check if the correct path is set.
For example, try the following:

```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

#### 6.1.2 If Fortran Compilation Fails

Please check if the `FC` variable in `src/homological_threading/fortran/Makefile` is set correctly.

### 6.2 Runtime Issues

#### 6.2.1 If Computation is Slow

- Increase the number of OpenMP threads: `export OMP_NUM_THREADS=8`
- Enable Python multiprocessing (`mp=True` option)

## 7. Developer Information

### 7.1 Extending the Code

To add new features, edit the following files:

- New analysis: `src/homological_threading/main.py`
- Modify input/output formats: `src/homological_threading/lammps_io.py`
- Speed up computation: `src/homological_threading/fortran/compute.f90`

### 7.2 Testing

To run tests:

```bash
python tests/test.py
