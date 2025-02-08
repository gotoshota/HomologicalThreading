import time
import argparse
import numpy as np
import MDAnalysis as mda
import h5py
import sys
import os

# モジュールのパスを追加
path_to_module = os.path.abspath(os.path.join(os.path.dirname(__file__), "../utils"))
sys.path.append(path_to_module)
import coordconv


def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="Parse LAMMPS data file path.")
    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Path to the INPUT LAMMPS data file, or tarfile.",
    )
    parser.add_argument(
        "--hdf", type=str, default="output.h5", help="Path to the OUTPUT HDF5 file."
    )
    parser.add_argument(
        "--lmpdata", type=str, default=None, help="Path to the OUTPUT data file."
    )
    return parser.parse_args()


def get_coords_by_mol(universe, mol):
    """
    Mol ID で指定された分子の座標を取得する
    """
    # filter atoms by mol id (same as resid id in MDAnalysis)
    atoms = universe.select_atoms(f"resid {mol}")
    coords = atoms.positions
    return coords


def replace_molids_in_data_file(original_data_file, output_data_file, N):
    """
    Replace molid in a LAMMPS data file. Each group of N atoms will have a unique molid.

    Parameters:
    original_data_file (str): Path to the original LAMMPS data file
    output_data_file (str): Path to the output data file with updated molids
    N (int): Number of atoms per molecule
    """
    with open(original_data_file, "r") as file:
        lines = file.readlines()

    with open(output_data_file, "w") as file:
        atom_section = False
        atom_id = 1
        for line in lines:
            if "Atoms" in line:  # Start of atom section
                atom_section = True
                file.write(line)
                continue

            # End of atom section if the line is empty or another section begins
            if atom_section and (
                "Bonds" in line or "Angles" in line or "Velocities" in line
            ):
                atom_section = False

            if atom_section and line.strip() != "":
                line_parts = line.split()
                molid = (atom_id - 1) // N + 1  # Calculate new molid
                line_parts[1] = str(molid)  # Replace molid in the line
                updated_line = " ".join(line_parts) + "\n"
                file.write(updated_line)
                atom_id += 1
            else:
                file.write(line)


if __name__ == "__main__":
    t1 = time.time()
    # Parse command line arguments
    args = parse_args()
    path_lmpdata = args.input
    path_h5data = args.hdf
    path_outdata = args.lmpdata
    # Check if the LAMMPS data file path is provided
    if not path_lmpdata:
        print(
            "Error: Please provide the path to the LAMMPS data file using the --input option."
        )
        exit()
    # Read LAMMPS Data file
    try:
        u = mda.Universe(path_lmpdata)
    except IOError:
        print("Error: Unable to read the LAMMPS data file. Please check the file path.")
        exit()
    # Get the number of beads, number of chains and total number of particles in system
    resids = u.atoms.resids
    nbeads = sum(1 for atom in u.atoms if atom.resid == 1)
    nchains = len(set(atom.resid for atom in u.atoms))
    nparticles = u.atoms.n_atoms
    if nbeads * nchains != nparticles:
        print("Error: The number of particles does NOT match the product of N*M")
        print(f"N    : the number of beads = {nbeads}")
        print(f"M    : the number of chains = {nchains}")
        exit()
    # WRAP
    box_dim = u.dimensions[0]  # Assume orthorhombic box
    center = np.array(
        [0.50 * box_dim, 0.50 * box_dim, 0.50 * box_dim]
    )  # the center of simulation cell
    # wrap the coordinates with conserve bond configuration
    wrapped_coords = coordconv.wrap_conserve_bondconf(
        center, box_dim, nparticles, nbeads, u.atoms.positions
    )
    a = coordconv.com(wrapped_coords)
    # put the center of mass of the system at the center of the simulation cell
    wrapped_coords = coordconv.transfer_min_dist(
        center, box_dim, nparticles, nbeads, wrapped_coords
    )
    # Update the coordinates in the Universe object
    u.atoms.positions = wrapped_coords
    a = coordconv.com(u.atoms.positions)
    # Save the coordinates to HDF5 file
    with h5py.File(path_h5data, "w") as f:
        f.create_dataset("/1_CoordConv/polywrap_coords", data=u.atoms.positions)
    # Save the coordinates to npy file
    # np.save(path_npydata, u.atoms.positions)
    # Write the updated universe to a LAMMPS data file
    if args.lmpdata is not None:
        with mda.Writer(path_outdata, u.atoms.n_atoms) as W:
            W.write(u)
        replace_molids_in_data_file(path_outdata, path_outdata, nbeads)
    t2 = time.time()
    print("Processing time: {:.2f} seconds".format(t2 - t1))
