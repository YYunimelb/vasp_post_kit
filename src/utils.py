import json
import os
import shutil
import numpy as np
from ase.io import read
import xml.etree.ElementTree as ET
import re







def structure_relaxation(path = ".",submit_command_std = "srun vasp > vaso"):
    shutil.copyfile("POSCAR","POSCAR-initial")
    # Start vasp relaxation
    os.system("rm -f track")
    with open("track","w") as f:
        f.writelines("Starting relaxation  \n")
        f.close()

    rel = "-1"

    while rel != "1":
        os.system(submit_command_std)
        with open("OSZICAR") as f:
            lines = f.readlines()
            last_line = lines[-1]
            rel = last_line.split()[0]

        with open("track","a") as f:
            f.writelines(rel+"\n")
            f.write(last_line)
        os.system("cp -f CONTCAR POSCAR")
        os.system(f"cp -f OUTCAR OUTCAR.{rel}")
        os.system(f"cp -f OSZICAR OSZICAR.{rel}")

    with open("OSZICAR") as f:
        lines = f.readlines()
        last_line = lines[-1]
        rel = last_line.split()[0]
    with open("track", "a") as f:
        f.writelines(rel)
        f.write(last_line)
    os.system('echo "Relaxation complete" >>track')

def mkdir(path):
    """
    :param path: str
        the path we need to make folder
    :return: none
    """

    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径\
def find_VBM_CBM_kpoints(path='OUTCAR'):
    '''

    This function are used to find the read the total number of electrons and the position of VBM and CBM

    Paremeters:

    path: str
        the path of OUTCAR,  # which should be the scf with high density of Kpoints

    return:
    NELECT:  int
        total number of valence electrons
    VBM_K[0][0]:  int
        the number k-points of VBM
    CBM_K[0][0]: int
        the number k-points of CBM

    for example
    :return  10 20  30
        10 valence electrons
        20th k points is the VBM
        30th k points is the CBM

    '''
    start = False
    find_k = False
    cord_x = []
    cord_y = []
    cord_z1 = []
    cord_z2 = []
    with open(path, 'r+') as f:
        for line in f.readlines():
            if line.find('vacuum') >= 0:
                Vline = float(line.split()[-1])
            if line.find('total number of electrons') >= 0:
                NELECT = int(float(line.split()[2]))  # NELECT  is the total number of electrons
            if (line.find('E-fermi') >= 0):
                start = True
            if start == True and (line.find('k-point') >= 0):
                cord = line.split()
                cord_x.append(float(cord[3]))
                cord_y.append(float(cord[4]))
                find_k = True
            if find_k == True and (line.split()[0] == f'{int(NELECT / 2)}'):
                cord = line.split()
                cord_z1.append(float(cord[1]))
            if find_k == True and (line.split()[0] == f'{int(NELECT / 2) + 1}'):
                cord = line.split()
                cord_z2.append((float(cord[1])))
                find_k = False
    f.close()
    VBM_K = np.where(cord_z1 == np.max(cord_z1))
    CBM_K = np.where(cord_z2 == np.min(cord_z2))
    return NELECT, int(VBM_K[0][0])+1, int(CBM_K[0][0])+1

def Modify_incar(file, label="IBRION", new_line="IBRION  = 2 "):
    # 读取INCAR文件
    with open(file, 'r') as f:
        lines = f.readlines()

    # 寻找并修改指定参数的值
    found = False
    for i in range(len(lines)):
        match = re.search(r'\b' + label + r'\b', lines[i])
        if match:
            lines[i] = new_line + '\n'
            found = True
            break

    # 如果未找到指定参数，则添加新行
    if not found:
        lines.append(new_line + '\n')

    # 写入修改后的内容
    with open(file, 'w') as f:
        f.writelines(lines)

def get_vaccum(LOCPOTfile = "LOCPOT",direction = "z"):
    from ase.calculators.vasp import VaspChargeDensity
    allowed = "xyzXYZ"
    if direction.islower():
        direction = direction.upper()

    vasp_charge = VaspChargeDensity(filename = LOCPOTfile)
    potl = vasp_charge.chg[-1]
    atoms = vasp_charge.atoms[-1]
    del vasp_charge

    if 'LOCPOT' in LOCPOTfile:
        potl=potl*atoms.get_volume()

    cell = atoms.cell

    # Find length of lattice vectors
    #--------------------------------
    latticelength = np.dot(cell, cell.T).diagonal()
    latticelength = latticelength**0.5

    # Read in potential data
    #------------------------
    ngridpts = np.array(potl.shape)
    totgridpts = ngridpts.prod()

    # Perform average
    #-----------------
    if direction=="X":
        idir = 0
        a = 1
        b = 2
    elif direction=="Y":
        a = 0
        idir = 1
        b = 2
    else:
        a = 0
        b = 1
        idir = 2
    a = (idir+1)%3
    b = (idir+2)%3
    # At each point, sum over other two indices
    average = np.zeros(ngridpts[idir],np.float64)
    for ipt in range(ngridpts[idir]):
        if direction=="X":
            average[ipt] = potl[ipt,:,:].sum()
        elif direction=="Y":
            average[ipt] = potl[:,ipt,:].sum()
        else:
            average[ipt] = potl[:,:,ipt].sum()

    if 'LOCPOT' in LOCPOTfile:
        # Scale by number of grid points in the plane.
        # The resulting unit will be eV.
        average /= ngridpts[a]*ngridpts[b]
    else:
        # Scale by size of area element in the plane,
        # gives unit e/Ang. I.e. integrating the resulting
        # CHG_dir file should give the total charge.
        area = np.linalg.det([ (cell[a,a], cell[a,b] ),
                               (cell[b,a], cell[b,b])])
        dA = area/(ngridpts[a]*ngridpts[b])
        average *= dA

    xdiff = latticelength[idir]/float(ngridpts[idir]-1)
    vaccum_energys = []
    for i in range(ngridpts[idir]):
        x = i*xdiff
        vaccum_energys.append(average[i])
    def get_average_round_the_max(vaccum_energys= [1,2,3,4,5,6,7,8]):
        max_index = vaccum_energys.index(max(vaccum_energys))
        sum = 0
        for i in range(5):
            sum += vaccum_energys[(max_index+i-2)%len(vaccum_energys)]
        return sum/5
    return get_average_round_the_max(vaccum_energys)
def collect_all_structures():
    """
    In the current directory, find all CONTCAR in subdirectory and collect the in "CONTCAR" folder
    :return:
    """
    mkdir("CONTCAR")
    for dirpath, dirnames, files in os.walk('.', topdown=False):
        for i in files:
            if i == "CONTCAR":
                a = dirpath.split("/")
                name = ""
                for i in a[1:]:
                    name = name + i+"_"
                shutil.copyfile(f"{dirpath}/CONTCAR",f"CONTCAR/{name}_CONTCAR")
def get_energy_from_OSZICAR(path = "OSZICAR"):
    """
    get the energy from OSZICAR
    :param path: str
        the path of OSZICAR
    :return: float
        the energy of this system
    """
    with open(path) as f:
        lines = f.readlines()
        a = lines[-1].split()
        return float(a[4])
def get_natoms_from_CONTCAR(path = "CONTCAR"):
    """
    Get the atoms number from CONTCAR
    :param path: srt
        the path of "CONTCAR"
    :return: int
        the number of atoms
    """
    with open(path) as f:
        lines = f.readlines()
        a = lines[6].split()
        natoms = 0
        for i in a:
            natoms += int(i)
        return  natoms
def collect_all_energys(path = "."):
    """
    In the current directory, find all OSZICAR in subdirectory and collect the in "CONTCAR" folder
    :return:
    """
    context = ""
    for dirpath, dirnames, files in os.walk('.', topdown=False):
        for i in files:
            if i == "OSZICAR":
                energy  = get_energy_from_OSZICAR(path = f"{dirpath}/OSZICAR")
                natoms = get_natoms_from_CONTCAR(path = f"{dirpath}/CONTCAR")
                name = ""
                a = dirpath.split("/")
                for i in a[1:]:
                    name = name + i.ljust(20)
                context = context + name +  str(energy).ljust(20)  + str(energy/natoms).ljust(20)    + "\n"
    with open("all_energy","w") as f:
        f.writelines(context)
def renew(path = "."):
    from fnmatch import fnmatchcase as match
    a = ["INCAR","POTCAR","KPOINTS","CONTCAR"]
    script = "*.sh"
    file_list = os.listdir(path)
    for i in file_list:
        if i not in a and not match(i,script):
            os.remove(f"{path}/{i}")
    os.system("cp CONTCAR POSCAR")
def check_track(path="track"):
    with open(path) as f:
        lines = f.readlines()
        a = lines[-1].split()
        if a[1] == "complete":
            b = lines[-2].split()
            if b[0] == "1":
                return True
        else:
            return False
def get_all_vasp_data_to_json(path ="."):
    import json
    total_information = []
    for dirpath, dirnames, files in os.walk(path, topdown=False):
        for i in files:
            if i == "OSZICAR":
                vasp = VaspData(dirpath)
                vasp.updata_all()
                total_information.append(json.dumps(vasp.__dict__))
    with open("Vasp.json","w") as f :
        json.dump(total_information,f)
def convert_to_excel():
    import pandas
def get_thickness(POSCAR_path="POSCAR"):
    def thickness_of_slab(atoms =[0.99,0.01]):
        def make_it_in_01(atoms=[2.554, -0.255]):
            new_atoms = []
            for atom in atoms:
                atom = atom % 1
                new_atoms.append(atom)
            return new_atoms

        def variance(atoms=[2.554, -0.255]):
            var = 0
            for atom in atoms:
                var += (atom - 0.5) ** 2
            return var / len(atoms)

        #structure = read(POSCAR_path)
        # atoms = structure.get_scaled_positions()[:,2]
        var = []
        for i in range(10):
            new_atoms = [j + i/10 for j in atoms]
            new_atoms = make_it_in_01(new_atoms)
            var.append(variance(new_atoms))
        index = var.index(min(var))
        new_atoms = [j + index / 10 for j in atoms]
        new_atoms = make_it_in_01(new_atoms)
        thickness = (max(new_atoms) - min(new_atoms))
        return thickness
    structure = read(POSCAR_path)
    vacuum = []
    for i in range(3):
        thick = thickness_of_slab(structure.get_scaled_positions()[:,i])
        vacuum.append(thick* structure.get_cell()[i,i])
    return vacuum
def read_eigenval(eigenval_path = "EIGENVAL"):
    kpoints = []
    weights = []
    energys = []
    with open(eigenval_path) as f:
        lines = f.readlines()
        for line in lines[7:]:
            line_split = line.split()
            if len(line_split) == 4:
                energy = []
                kpoints.append([float(line_split[0]),float(line_split[1]),float(line_split[2])])
                weights.append(float(line_split[3]))
            elif len(line_split) == 3:
                energy.append((float(line_split[1]),float(line_split[2])))
            else:
                energys.append(energy)
        energys.append(energy)
    return kpoints,weights,energys
def cal_rec_vector(lattvec):
    def DotProduct(v1, v2):
        dotproduct = 0.0
        for i in range(0, 3):
            dotproduct = dotproduct + v1[i] * v2[i]
        return dotproduct

    def CrossProduct(v1, v2):
        # v3 = v1 WILL LEAD TO WRONG RESULT
        v3 = []
        v3.append(v1[1] * v2[2] - v1[2] * v2[1])
        v3.append(v1[2] * v2[0] - v1[0] * v2[2])
        v3.append(v1[0] * v2[1] - v1[1] * v2[0])
        return v3

    def CalcRecVectors(lattvec):
        pi = 3.141592653589793
        a1 = lattvec[0]
        a2 = lattvec[1]
        a3 = lattvec[2]
        b1 = CrossProduct(a2, a3)
        b2 = CrossProduct(a3, a1)
        b3 = CrossProduct(a1, a2)
        volume = DotProduct(a1, CrossProduct(a2, a3))
        RecVec = [b1, b2, b3]
        # it follows the definition for b_j: b_i * b_j = 2pi * delta(i,j)
        for i in range(0, 3):
            for j in range(0, 3):
                #  RecVec[i][j] = RecVec[i][j] * 2 * pi / volume   this one is true but different from vasp case
                RecVec[i][j] = RecVec[i][j]  / volume  # The result of this one is consistent with VASP
        return RecVec
    return CalcRecVectors(lattvec)

