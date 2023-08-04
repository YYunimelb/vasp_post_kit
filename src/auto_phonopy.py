import os
import shutil
import vasp.utils as utils
from vasp.environment import All_Environment_Parameter

# 设置Vasp执行命令
all_parameter = All_Environment_Parameter()
submit_command_std = all_parameter.submit_vasp_command_std



class Phonon():
    def __init__(self,config = {} ,envir_config = {submit_command_std:" "}):
        self.pwd = os.getcwd()
        self.structure_path = os.path.join(self.pwd, "structure")
        self.phonopy_path = os.path.join(self.pwd, "phonopy")

    # 检查文件是否存在
    def __file_exists(self, file_path):
        return os.path.exists(file_path)

    # 计算晶格常数
    def __calculate_the_lattice_constant(self, a):
        return ((float(a[0])) ** 2 + (float(a[1])) ** 2 + (float(a[2]) ** 2)) ** 0.5

    # 从POSCAR文件获取晶格常数
    def __get_lattice_constant_from_poscar(self, path_of_POSCAR="POSCAR"):
        if self.__file_exists(path_of_POSCAR):
            with open(path_of_POSCAR) as f:
                lines = f.readlines()
                a = lines[2].split()
                b = lines[3].split()
                c = lines[4].split()
                return (self.__calculate_the_lattice_constant(a),
                        self.__calculate_the_lattice_constant(b),
                        self.__calculate_the_lattice_constant(c))

    # 获取超晶胞大小
    def get_supercell_size(self, supercell_lattice=[15, 15, 15]):
        a, b, c = self.__get_lattice_constant_from_poscar(os.path.join(self.structure_path, "CONTCAR"))
        self.supercell = [int(supercell_lattice[0] // a + 1),
                          int(supercell_lattice[1] // b + 1),
                          int(supercell_lattice[2] // c + 1)]
        print(self.supercell)

    # 生成超晶胞
    def make_supercell(self):
        utils.mkdir(self.phonopy_path)
        shutil.copyfile(os.path.join(self.structure_path, "CONTCAR"), os.path.join(self.phonopy_path, "POSCAR"))
        shutil.copyfile(os.path.join(self.structure_path, "POTCAR"), os.path.join(self.phonopy_path, "POTCAR"))
        shutil.copyfile(os.path.join(self.structure_path, "INCAR"), os.path.join(self.phonopy_path, "INCAR"))
        shutil.copyfile(os.path.join(self.structure_path, "INCAR"), os.path.join(self.phonopy_path, "INCAR1"))
        os.chdir(self.phonopy_path)
        os.system(f'phonopy -d --dim="{self.supercell[0]} {self.supercell[1]} {self.supercell[2]}"')

    # 生成K点
    def make_kpoints(self, kpoints=[3, 3, 1]):
        content = f"""
        pymatgen v2020.8.13 with grid density = 654 / number of atoms
        0
        Gamma
        {kpoints[0]} {kpoints[1]} {kpoints[2]}
        """
        with open(os.path.join(self.phonopy_path, "KPOINTS"), "w") as f:
            f.writelines(content)

    # 创建band.conf文件
    def make_band_conf(self):
        content = f"""
        ATOM_NAME =haha
        DIM = {self.supercell[0]} {self.supercell[1]} {self.supercell[2]}
        PRIMITIVE_AXES=Auto
        BAND =0.0 0.0 0.0  0.5 0.0 0.0  0.5  0.5   0.0   0.0  0.0  0.0  0.0  0.5  0.0
        BAND_POINTS = 505"""
        with open(os.path.join(self.phonopy_path, "band.conf"), "w") as f:
            f.writelines(content)

    # 修改INCAR文件
    def modify_incar(self):
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "ISIF", " ISIF = 2")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "NSW", " NSW = 1")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "IBRION", " IBRION = -1")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "EDIFFG", "# EDIFFG  =  -0.001")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "LCHARG", " LCHARG = .TRUE. ")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "LWAVE", " LWAVE = .TRUE. ")
        utils.Modify_incar(os.path.join(self.phonopy_path, "INCAR"), "EDIFF", "  EDIFF  =  1E-08 ")

    # 执行001任务
    def do_001(self):
        utils.mkdir(os.path.join(self.phonopy_path, "001"))
        shutil.copyfile(os.path.join(self.phonopy_path, "INCAR"), os.path.join(self.phonopy_path, "001", "INCAR"))
        shutil.copyfile(os.path.join(self.phonopy_path, "POTCAR"), os.path.join(self.phonopy_path, "001", "POTCAR"))
        shutil.copyfile(os.path.join(self.phdo_all()
