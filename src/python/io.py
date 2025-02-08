class lammpsData:
    """
    A class to read and store LAMMPS data from a file.

    Attributes:
        filename (str): The name of the file containing LAMMPS data.
        atoms (list): A list to store atom data.
        bonds (list): A list to store bond data.
        angles (list): A list to store angle data.
        dihedrals (list): A list to store dihedral data.
        impropers (list): A list to store improper data.
        masses (list): A list to store mass data.
        box (dict): A dictionary to store box dimensions.
    """

    def __init__(self, filename):
        """
        Initializes the lammpsData object by reading data from the specified file.

        Args:
            filename (str): The name of the file containing LAMMPS data.
        """
        self.filename = filename
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.masses = []
        self.atom_types = 0
        self.bond_types = 0
        self.angle_types = 0
        self.dihedral_types = 0
        self.improper_types = 0
        self.box = {}
        self.readData()

    def readData(self):
        """
        Reads data from the file and stores it in the appropriate attributes.
        """
        with open(self.filename, "r") as f:
            lines = [line.strip() for line in f.readlines()]  # 余分な空白や改行を除去
            section_map = {
                "Atoms": self.atoms,
                "Bonds": self.bonds,
                "Angles": self.angles,
                "Dihedrals": self.dihedrals,
                "Impropers": self.impropers,
                "Masses": self.masses,
            }
            type_map = {
                "atom types": "atom_types",
                "bond types": "bond_types",
                "angle types": "angle_types",
                "dihedral types": "dihedral_types",
                "improper types": "improper_types",
            }

            for i, line in enumerate(lines):
                for type_name, attr_name in type_map.items():
                    if type_name in line:
                        setattr(
                            self, attr_name, int(line.split()[0])
                        )  # クラス変数を直接更新

                for section, storage in section_map.items():
                    if line.startswith(section):
                        j = i + 2  # セクションヘッダーの次の行からデータを取得
                        while j < len(lines):
                            current_line = lines[j].strip()
                            if not current_line or current_line.startswith(
                                "#"
                            ):  # 空行やコメントをスキップ
                                j += 1
                                continue
                            if any(
                                current_line.startswith(next_section)
                                for next_section in section_map.keys()
                            ):
                                break
                            storage.append(current_line.split())  # 適切なデータのみ追加
                            j += 1

                if "xlo" in line:
                    self.box["x"] = (float(line.split()[0]), float(line.split()[1]))
                if "ylo" in line:
                    self.box["y"] = (float(line.split()[0]), float(line.split()[1]))
                if "zlo" in line:
                    self.box["z"] = (float(line.split()[0]), float(line.split()[1]))

    def writeData(self, filename):
        """
        Writes data to a file.

        Args:
            filename (str): The name of the file to write data to.
        """
        with open(filename, "w") as f:
            f.write("LAMMPS data by LammpsDumpReader\n\n")
            f.write(f"{len(self.atoms)} atoms\n")
            f.write(f"{len(self.bonds)} bonds\n")
            f.write(f"{len(self.angles)} angles\n")
            f.write(f"{len(self.dihedrals)} dihedrals\n")
            f.write(f"{len(self.impropers)} impropers\n\n")
            f.write(f"{self.atom_types} atom types\n")
            f.write(f"{self.bond_types} bond types\n")
            f.write(f"{self.angle_types} angle types\n")
            f.write(f"{self.dihedral_types} dihedral types\n")
            f.write(f"{self.improper_types} improper types\n\n")
            f.write(f"{self.box['x'][0]} {self.box['x'][1]} xlo xhi\n")
            f.write(f"{self.box['y'][0]} {self.box['y'][1]} ylo yhi\n")
            f.write(f"{self.box['z'][0]} {self.box['z'][1]} zlo zhi\n\n")
            f.write("Masses\n\n")
            for mass in self.masses:
                f.write(" ".join(mass) + "\n")
            f.write("\nAtoms\n\n")
            for atom in self.atoms:
                f.write(" ".join(atom) + "\n")
            f.write("\nBonds\n\n")
            for bond in self.bonds:
                f.write(" ".join(bond) + "\n")
            f.write("\nAngles\n\n")
            for angle in self.angles:
                f.write(" ".join(angle) + "\n")
            f.write("\nDihedrals\n\n")
            for dihedral in self.dihedrals:
                f.write(" ".join(dihedral) + "\n")
            f.write("\nImpropers\n\n")
            for improper in self.impropers:
                f.write(" ".join(improper) + "\n")

    def polymerWrap(self):
        """
        Wraps coordinates of polymer chains in the box
        to conserve bond configuration.
        """
        for i, atom in enumerate(self.atoms):
            if i == 0:
                molecule_id = int(atom[0])
            elif atom[0] == molecule_id:
                atom[2] = str(float(atom[2]) % (self.box["x"][1] - self.box["x"][0]))
                atom[3] = str(float(atom[3]) % (self.box["y"][1] - self.box["y"][0]))
                atom[4] = str(float(atom[4]) % (self.box["z"][1] - self.box["z"][0]))
            else:
                molecule_id += 1

    def __str__(self):
        """
        Returns a string representation of the lammpsData object.

        Returns:
            str: A string representation of the lammpsData object.
        """
        return f"lammpsData({self.filename})"

    def __repr__(self):
        """
        Returns a string representation of the lammpsData object.

        Returns:
            str: A string representation of the lammpsData object.
        """
        return f"lammpsData({self.filename})"


if __name__ == "__main__":
    data = lammpsData("../../../murashimaPolym/N2000/msd.N2000.1.data")
    print(f"{len(data.atoms)=}")
    print(f"{len(data.bonds)=}")
    print(f"{data.bonds[0]=}")
    print(f"{len(data.angles)=}")
    print(f"{len(data.dihedrals)=}")
    print(f"{len(data.impropers)=}")
    print(f"{len(data.masses)=}")
    print(f"{data.atom_types=}")
    print(f"{data.bond_types=}")
    print(f"{data.angle_types=}")
    print(f"{data.dihedral_types=}")
    print(f"{data.improper_types=}")
    print(f"{data.box=}")
    print(f"{data.box['x']=}")
    print(f"{data.box['x'][0]=}")
    a = data.atoms[0][2]
    data.polymerWrap()
    print(f"{data.atoms[0][2] - a}")
    data.writeData("test.data")
