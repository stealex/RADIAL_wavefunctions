import configparser
import subprocess
import numpy as np
import os
import shutil
from mendeleev import element
import re
from python_util import dhfs_util
import sys


def main():
    base_dir = os.getcwd()
    print("base dir = ", base_dir)

    final_dir = "binding_energies"
    try:
        os.mkdir(final_dir)
    except OSError as error:
        print(error)

    shutil.copy2(base_dir+"/build/DHFS", final_dir+"/DHFS")
    shutil.copy2(base_dir+"/build/constants.mod", final_dir+"/constants.mod")

    dhfs_input_repo = base_dir+"/DHFS/input"

    os.chdir(final_dir)
    be_file = open("binding_energies.dat", 'w')
    be_file.write("Z,HolePosition,BindingEnergy,GSBindingEnergy,Delta\n")
    for iz in range(23, 24):
        gs_binding_energy = 1.E9
        # run dhfs with gs config of the final atom
        shutil.copy2(dhfs_input_repo+"/"+f'Z{iz:0>3d}.in', "tmp.in")
        tmp_file = open("tmp.in", 'r')
        subprocess.run("./DHFS", stdin= tmp_file)
        tmp_file.close()
        
        # extract gs energy of final atom
        out_file = open("dhfs.dat", 'r')
        lines = out_file.readlines()
            
        for line in lines:
            tokens = line.split()
            if (len(tokens) < 2):
                continue
            if (tokens[0] == "Total") & (tokens[1] == "binding"):
                gs_binding_energy = float(tokens[7])
        

        # prepare run with holes

        # store header and footer lines 
        # number of shells has to be taken from the initial atom
        file_header = open(dhfs_input_repo+"/"+f'Z{iz:0>3d}.in', 'r')
        lines = file_header.readlines()
        header_lines = []
        footer_lines = []
        for line in lines:
            tokens = line.split()
            if (tokens[0] == "C1") | (tokens[0] == "C2"):
                header_lines.append(line)
            if (tokens[0] == "C4") | (tokens[0] == "C5") | (tokens[0] == "C6") | (tokens[0] == "C7"):
                footer_lines.append(line)
        file_header.close()
        # print("header lines = ", header_lines)
        # print("footer lines = ", footer_lines)

        # store initial electron configuration
        file_initial = open(dhfs_input_repo+"/"+f'Z{iz+1:0>3d}.in', 'r')
        lines = file_initial.readlines()
        electron_config_lines = []
        n_shells = ""
        for line in lines:
            tokens = line.split()
            if (tokens[0] == "C3"):
                electron_config_lines.append(line)
            if (tokens[0] == "C2"):
                tmp = tokens[-9]
                n_shells = tmp.split(")")[0]
        file_initial.close()
        # print("nShells = ", n_shells)
        # print("electron config = ", electron_config_lines)
        # replace n_shells in header

        c2_line = re.sub(r' \d{1,2}\)', n_shells+")", header_lines[1])
        print("c2_line = ", c2_line)
        header_lines[1] = c2_line
        
        # write new file with hole in K, L1, L2, L3 shells
        hole_pos_dict={0: "K", 1: "L1", 2: "L2", 3: "L3"}
        for iHole in range(4):
            file_final = open(f'Z{iz:0>3d}_{hole_pos_dict[iHole]:s}hole.in', 'w')

            # write header lines
            for line in header_lines:
                file_final.write(line)
            
            for iShell in range(len(electron_config_lines)):
                line = electron_config_lines[iShell]
                # modify line: put a hole in it
                if iShell == iHole:
                    tokens = line.split()  
                    occ_num = float(tokens[5].split(")")[0])
                    tokens[5] = f'{occ_num-1:6.3f}'
                    new_line = f'C3 ( {tokens[2]} {tokens[3]} {tokens[4]}{tokens[5]})\n'
                    file_final.write(new_line)
                else:
                    file_final.write(line)
            
            # write footer lines
            for line in footer_lines:
                file_final.write(line)
            file_final.close()

            # run dhfs
            file_final = open(f'Z{iz:0>3d}_{hole_pos_dict[iHole]:s}hole.in', 'r')
            subprocess.run("./DHFS", stdin= file_final)
            file_final.close()

            # extract result
            out_file = open("dhfs.dat", 'r')
            lines = out_file.readlines()
            for line in lines:
                tokens = line.split()
                if (len(tokens) < 2):
                    continue
                if (tokens[0] == "Total") & (tokens[1] == "binding"):
                    binding_energy = float(tokens[7])
                    tmp = f'{iz:d},{iHole:d},{binding_energy:14.7e},{gs_binding_energy:14.7e},{binding_energy-gs_binding_energy:14.7e}\n'
                    print("fucking tmp = ", tmp)
                    be_file.write(tmp)

            out_file.close()
    be_file.close()
    
if __name__ == "__main__":
    main()