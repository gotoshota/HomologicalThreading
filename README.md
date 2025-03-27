# Homological Threading Documentation

- [Japanese Version](README.ja.md)
- [English Version](README.md)

# Homological Threading

This project is designed to quantify the threading (entanglement) of cyclic polymers using persistent homology. By applying persistent homology, the program analyzes the geometric and topological features (e.g., holes, voids, connected components) in a multi-scale manner to capture the complex intertwinement between cyclic polymers.

---

## Table of Contents
- [Homological Threading Documentation](#homological-threading-documentation)
- [Homological Threading](#homological-threading)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Directory Structure](#directory-structure)
  - [3. System Requirements and Dependencies](#3-system-requirements-and-dependencies)
    - [3.1 Requirements](#31-requirements)
      - [3.1.1 Python Environment](#311-python-environment)
      - [3.1.2 Fortran Compiler](#312-fortran-compiler)
      - [3.1.3 CGAL (Computational Geometry Algorithms Library)](#313-cgal-computational-geometry-algorithms-library)
  - [4. Installation Instructions](#4-installation-instructions)
    - [4.1 Installing Dependencies](#41-installing-dependencies)
    - [4.2 Setting Up the Python Virtual Environment](#42-setting-up-the-python-virtual-environment)
    - [4.3 Building the Project](#43-building-the-project)
  - [5. Usage](#5-usage)
    - [5.1 Calculation and Analysis of Persistence Diagrams](#51-calculation-and-analysis-of-persistence-diagrams)
    - [5.2 Betti Numbers Calculation](#52-betti-numbers-calculation)
    - [5.3 Visualization of Results](#53-visualization-of-results)
  - [6. Running Tests](#6-running-tests)
  - [7. Troubleshooting](#7-troubleshooting)
    - [7.1 Installation Issues](#71-installation-issues)
    - [7.2 Runtime Issues](#72-runtime-issues)
  - [8. Developer Information](#8-developer-information)
    - [8.1 Extending the Code](#81-extending-the-code)
    - [8.2 Changes to Repository Structure](#82-changes-to-repository-structure)
  - [9. License](#9-license)

---

## Features
- **Persistent Homology Analysis**: Analyzes point cloud data through filtrations at multiple scales to quantify topological features such as connected components, holes, and voids.
- **Threading Quantification**: Evaluates the entanglement between single and multiple cyclic polymers using persistence diagrams and Betti numbers.
- **High-Performance Computation**: The computational core is implemented in Fortran, ensuring efficient processing even for large-scale datasets.

---

## Directory Structure

The main files and directories of the project are as follows:

```
HomologicalThreading/          ← Project root
├── .gitignore
├── .python-version
├── build.sh                   ← Build script
├── main.py                    ← Entry point (interface)
├── pyproject.toml             ← Project configuration file
├── README.md                  ← This English README
├── README.ja.md               ← Japanese README
├── uv.lock                    ← Dependency lock file
├── data/                     ← Sample data, etc.
│   ├── N10M100.data          ← Sample LAMMPS data file
│   └── (other data)
├── docs/                     ← Documentation
│   ├── build/
│   └── source/
│       ├── _static/
│       └── _templates/
├── scripts/                  ← Analysis and visualization scripts
│   ├── __init__.py
│   ├── analysis.py           ← Persistence diagram and Betti numbers analysis
│   └── plot.py               ← Visualization (integrated plotting functions)
├── src/                      ← Source code
│   └── homological_threading/
│       ├── __init__.py
│       ├── lammps_io.py      ← Input/output for LAMMPS data files
│       ├── main.py           ← Core functionality implementation
│       └── fortran/          ← Fortran implementation for high-performance computation
│           ├── __init__.py
│           ├── compute.f90   ← Core routines for threading and Betti numbers calculation
│           └── Makefile
└── tests/                    ← Test code
    └── test.py
```

*Note: The directory structure is based on the latest project layout.*

---

## 3. System Requirements and Dependencies

### 3.1 Requirements

#### 3.1.1 Python Environment
- Python 3.13 or higher
- uv (Python project management tool)

Install uv using:
```bash
pip install uv
```

#### 3.1.2 Fortran Compiler
- gfortran (GNU Fortran Compiler) or ifx (Intel Fortran Compiler)

On macOS:
```bash
brew install gcc
```

#### 3.1.3 CGAL (Computational Geometry Algorithms Library)
CGAL is required as a dependency for the HomCloud library.

On macOS:
```bash
brew install cgal
```

Alternatively, compile BOOST and CGAL from source for a local installation.

---

## 4. Installation Instructions

### 4.1 Installing Dependencies
Ensure that the above dependencies (Python, Fortran compiler, CGAL, etc.) are installed.

### 4.2 Setting Up the Python Virtual Environment

Use uv to create and synchronize the virtual environment:
```bash
uv sync
```

If the CGAL include path is required, specify it as follows:
```bash
CPLUS_INCLUDE_PATH=/path/to/CGAL-5.6.2/include:$CPLUS_INCLUDE_PATH uv sync
```

### 4.3 Building the Project

Build the Fortran code and perform any other necessary build steps with:
```bash
./build.sh
```

If the Fortran compiler is not detected, adjust the `FC` variable in `src/homological_threading/fortran/Makefile` appropriately.

---

## 5. Usage

### 5.1 Calculation and Analysis of Persistence Diagrams
Compute the persistence diagrams (PD) of cyclic polymers from a LAMMPS data file and analyze threading:
```bash
uv run python main.py analysis pd -i data/N10M100.data -o output_directory
```
This command:
1. Reads cyclic polymer coordinates from the LAMMPS data file.
2. Computes the persistence diagram (PD_i) for individual polymers.
3. Analyzes threading (e.g., by computing PD_i_cup_j) between multiple polymers.
4. Saves the analysis results in HDF5 format.

### 5.2 Betti Numbers Calculation

Calculate Betti numbers from the saved HDF5 files:
```bash
uv run python main.py analysis betti -i output_directory/*.h5 -f output_directory/analysis.h5
```

### 5.3 Visualization of Results

Visualize the analysis results using:
```bash
uv run python main.py plot -i output_directory/analysis.h5
```

---

## 6. Running Tests

Run the test suite with:
```bash
python tests/test.py
```

---

## 7. Troubleshooting

### 7.1 Installation Issues

- **HomCloud Build Failure**: Ensure that the CGAL include path is set correctly.
- **Fortran Compilation Failure**: Check that the `FC` variable in `src/homological_threading/fortran/Makefile` is correctly configured.

### 7.2 Runtime Issues

- **Slow Computation**:
  - Increase the number of threads by setting the `OMP_NUM_THREADS` environment variable (e.g., `export OMP_NUM_THREADS=8`).
  - Consider enabling Python multiprocessing.

---

## 8. Developer Information

### 8.1 Extending the Code
To add new features, edit:
- Core analysis logic: `src/homological_threading/main.py`
- Input/output routines: `src/homological_threading/lammps_io.py`
- High-performance routines: `src/homological_threading/fortran/compute.f90` and `Makefile`
- Tests: `tests/test.py`

### 8.2 Changes to Repository Structure
Update root-level files (e.g., `main.py`, `build.sh`, `pyproject.toml`) as needed.

---

## 9. License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
