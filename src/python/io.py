class LammpsData:
    """
    LAMMPSのデータファイルを読み込み，データを格納するクラス

    Attributes:
        filename (str): LAMMPSデータファイルのパス
        box (Box): ボックス情報を保持するオブジェクト
        atoms (Atom): 原子情報を保持するオブジェクト
        bonds (Bond): ボンド情報を保持するオブジェクト
        angles (Angle): 角度情報を保持するオブジェクト
        現在，2面角 (Dihedrals) や不整合角 (Impropers) は未対応． また，将来の拡張も考慮していない．
    """

    def __init__(self, filename=None):
        """
        LammpsDataオブジェクトを初期化し，ファイルが指定されていればデータを読み込む

        Args:
            filename (str): LAMMPSデータファイルのパス
        """
        self.filename = filename
        self.atoms = self.Atom()  # 一貫して atoms としています
        self.bonds = self.Bond()
        self.angles = self.Angle()
        self.box = self.Box()

        if filename is not None:
            self.read(filename)

    class Box:
        """
        ボックス情報を保持するクラス

        Attributes:
            x (tuple): x軸の範囲 (xlo, xhi)
            y (tuple): y軸の範囲 (ylo, yhi)
            z (tuple): z軸の範囲 (zlo, zhi)
            lx (float): x軸の長さ
            ly (float): y軸の長さ
            lz (float): z軸の長さ
        """

        def __init__(self):
            self.x = (0.0, 0.0)
            self.y = (0.0, 0.0)
            self.z = (0.0, 0.0)
            self.lx = 0.0
            self.ly = 0.0
            self.lz = 0.0

    class Atom:
        """
        原子情報を保持するクラス

        Attributes:
            id (list of int): 原子IDのリスト
            mol_id (list of int): 分子IDのリスト
            type (list of int): 原子タイプのリスト
            coords (list of tuple): 座標 (x, y, z) のリスト
            num_atoms (int): 原子の総数
        """

        def __init__(self):
            self.id = []
            self.mol_id = []
            self.type = []
            self.coords = []
            self.num_atoms = 0
            self.num_mols = 0
            self.num_types = 0

    class Bond:
        """
        ボンド情報を保持するクラス

        Attributes:
            id (list of int): ボンドIDのリスト
            type (list of int): ボンドタイプのリスト
            atoms (list of tuple): ボンドを構成する原子のリスト
            num_bonds (int): ボンドの総数
            num_types (int): ボンドタイプの数
        """

        def __init__(self):
            self.id = []
            self.type = []
            self.atoms = []
            self.num_bonds = 0
            self.num_types = 0

    class Angle:
        """
        角度情報を保持するクラス

        Attributes:
            id (list of int): 角度IDのリスト
            type (list of int): 角度タイプのリスト
            atoms (list of tuple): 角度を構成する原子のリスト
            num_angles (int): 角度の総数
            num_types (int): 角度タイプの数
        """

        def __init__(self):
            self.id = []
            self.type = []
            self.atoms = []
            self.num_angles = 0
            self.num_types = 0

    def read(self, filename):
        """
        LAMMPSデータファイルからデータを読み込み，各属性に格納する

        Args:
            filename (str): 読み込むファイル名
        """
        if filename is not None:
            self.filename = filename

        with open(self.filename, "r") as f:
            lines = f.readlines()

        # --- 1. ヘッダー部からボックス情報や型数を取得 ---
        for line in lines:
            line_strip = line.strip()
            if "xlo" in line_strip and "xhi" in line_strip:
                parts = line_strip.split()
                self.box.x = (float(parts[0]), float(parts[1]))
                self.box.lx = float(parts[1]) - float(parts[0])
            elif "ylo" in line_strip and "yhi" in line_strip:
                parts = line_strip.split()
                self.box.y = (float(parts[0]), float(parts[1]))
                self.box.ly = float(parts[1]) - float(parts[0])
            elif "zlo" in line_strip and "zhi" in line_strip:
                parts = line_strip.split()
                self.box.z = (float(parts[0]), float(parts[1]))
                self.box.lz = float(parts[1]) - float(parts[0])
            elif "atoms" in line_strip:
                parts = line_strip.split()
                self.atoms.num_atoms = int(parts[0])
            elif "bonds" in line_strip:
                parts = line_strip.split()
                self.bonds.num_bonds = int(parts[0])
            elif "angles" in line_strip:
                parts = line_strip.split()
                self.angles.num_angles = int(parts[0])
            elif "atom types" in line_strip:
                parts = line_strip.split()
                self.atoms.num_types = int(parts[0])
            elif "bond types" in line_strip:
                parts = line_strip.split()
                self.bonds.num_types = int(parts[0])
            elif "angle types" in line_strip:
                parts = line_strip.split()
                self.angles.num_types = int(parts[0])

        # --- 2. セクション毎にデータをパース ---
        section_names = ["Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"]
        current_section = None
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # セクションヘッダーの検出
            if any(line.startswith(sec) for sec in section_names):
                # セクション名を現在のセクションとして記憶
                current_section = line.split()[0]
                # セクション名の行とその次の空行（またはコメント行）をスキップ
                i += 2
                continue

            # セクション外なら次の行へ
            if current_section is None:
                i += 1
                continue

            # 空行やコメント行はスキップ
            if not line or line.startswith("#"):
                i += 1
                continue

            # 新たなセクションが始まった場合は current_section をリセット
            if any(line.startswith(sec) for sec in section_names):
                current_section = None
                continue

            # データ行のパース
            parts = line.split()
            if current_section == "Atoms":
                # 例: "1 1 1 0.0 0.0 0.0 ..." (原子ID, 分子ID, タイプ, x, y, z, ...)
                if len(parts) >= 6:
                    self.atoms.id.append(int(parts[0]))
                    self.atoms.mol_id.append(int(parts[1]))
                    self.atoms.type.append(int(parts[2]))
                    self.atoms.coords.append(
                        (float(parts[3]), float(parts[4]), float(parts[5]))
                    )
            # 他のセクション (Masses, Bonds, など) のパース処理も同様に追加可能

            i += 1
        # 被ってない，mol_id の数を数える
        self.atoms.num_mols = len(set(self.atoms.mol_id))

    def __str__(self):
        return f"LammpsData({self.filename})"

    def __repr__(self):
        return self.__str__()


# 動作確認用（必要に応じてパスを適宜変更してください）
if __name__ == "__main__":
    data = LammpsData("../../../murashimaPolym/N2000/msd.N2000.1.data")
    # Test box
    print(f"{data.box.x=}")
    print(f"{data.box.y=}")
    print(f"{data.box.z=}")
    print(f"{data.box.lx=}")
    print(f"{data.box.ly=}")
    print(f"{data.box.lz=}")

    # Test atoms
    print(f"{data.atoms.num_atoms=}")
    print(f"{data.atoms.num_mols=}")
    print(f"{data.atoms.num_types=}")

    # Test bonds
    print(f"{data.bonds.num_bonds=}")
    print(f"{data.bonds.num_types=}")
