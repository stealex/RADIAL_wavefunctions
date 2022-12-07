import configparser
import subprocess
import numpy as np
import os
import shutil
from mendeleev import element
import re

def move_dhfs_output(z:int, nuc):
	if (nuc != "initial") and ("final" not in nuc):
		print("nuc parameter can be initial or contain daugher")
		return
	os.system(f'mv dhfs{z:0>3}.tab dhfs_{nuc:s}.tab')
	os.system(f'mv dhfs.dat dhfs_{nuc:s}.dat')
	os.system(f'mv potden.dat potden_{nuc:s}.dat')
	os.system(f'mv rwavefcts.dat rwavefcts_{nuc:s}.dat')

def make_potential_from_dhfs(nuc):
	# read potential for initial nucleus
	with open(f"potden_{nuc}.dat") as potden_file:
		radius_tmp = []
		rv_nuc_tmp = []
		rv_el_tmp = []
		rv_ex_tmp = []
		rho_tmp = []
		for line in potden_file.readlines():
			if "#" in line:
				continue
			tokens = line.strip().split(" ")
			tokens = list(filter(lambda a: a != '', tokens))
			radius_tmp.append(float(tokens[0]))
			rv_nuc_tmp.append(float(tokens[2]))
			rv_el_tmp.append(float(tokens[3]))
			rv_ex_tmp.append(float(tokens[4]))
			rho_tmp.append(float(tokens[5]))
		
		radius = np.array(radius_tmp)
		rv_nuc = np.array(rv_nuc_tmp)
		rv_el = np.array(rv_el_tmp)
		rv_ex = np.array(rv_ex_tmp)
		rho = np.array(rho_tmp)

		potential_with_Latter = rv_nuc+rv_el+rv_ex
		potential_no_Latter = rv_nuc+rv_el - 1.5*((3./np.pi)**(1./3.))*radius*(rho**(1./3.))
		
		np.savetxt(f'potential_with_Latter_{nuc}.dat', np.c_[radius, potential_with_Latter])
		np.savetxt(f'potential_no_Latter_{nuc}.dat', np.c_[radius, potential_no_Latter])


project_base_dir = os.getcwd()

config = configparser.ConfigParser()
config.read("dhfs.ini")

# initial nucleus name and the process it undergoes
initial_nucleus = config["DIRAC"]["nucleusName"]
process = config["DIRAC"]["processName"]
z_initial = int(config["DIRAC"]["zInitial"])
initial_element = element(z_initial)

# create directory to store results
final_dir = initial_nucleus+"_"+process
try:
	os.mkdir(final_dir)
except OSError as error:
	print(error)

shutil.copy2(project_base_dir+"/build/DHFS", final_dir+"/DHFS")
shutil.copy2(project_base_dir+"/build/constants.mod", final_dir+"/constants.mod")


# copy the dhfs .in file to the final directory
if int(config["DHFS"]["USE_STANDARD"]) == 1:
	repo = project_base_dir+"/"+config["DHFS"]["REPO"]
	shutil.copyfile(repo+"/"+f'Z{z_initial:0>3d}.in', final_dir+"/dhfs_initial.in")
else:
	shutil.copyfile(str(config["INPUT_FILE"]), final_dir+"/dhfs_initial.in")

os.chdir(final_dir)

# Run DHFS with input file specified in config and rename output files
dhfs_input_file = open("dhfs_initial.in")
subprocess.run("./DHFS", stdin= dhfs_input_file)
dhfs_input_file.close()
dhfs_input_file = open("dhfs_initial.in")
move_dhfs_output(z_initial, "initial")

# Create dhfs input files for final nucleus based on process and run DHFS
if (process == "EC") or (process == "2EC"):
	print("EC")
elif (process == "betaMinus") or (process == "2betaMinus"):
	print(process)
	# only 1 final configuration
	z_final = z_initial+1
	if (process == "2betaMinus"):
		z_final = z_initial+2

	final_element = element(z_final)
	dhfs_initial_lines = dhfs_input_file.readlines()
	print("Writting dhfs input file for final nulcues ... ")
	with open("dhfs_final.in", "w") as dhfs_final_file:
		for line in dhfs_initial_lines:
			if "C1" in line or "C2" in line:
				line = line.replace(f'{z_initial: >3d}', f'{z_final: >3d}')
				line = line.replace("Neutral atom", "Positive ion")
				line = line.replace(str(initial_element.symbol), str(final_element.symbol))
				line = line.replace(str(initial_element.name), str(final_element.name))
			
			if "C4" in line:
				initial_weight = re.findall(r"\d+\.\d+", line)
				print("initial weight = ", initial_weight)
				final_weight = final_element.mass
				print("final weight = ", final_weight)
				line = line.replace(f'{float(initial_weight[0]): >9.4f}', f'{float(final_weight): >9.4f}')
				print(line)
			dhfs_final_file.write(line)

	dhfs_final_file = open("dhfs_final.in")
	subprocess.run("./DHFS", stdin= dhfs_final_file)
	dhfs_final_file.close()

	move_dhfs_output(z_final, "final")
	print("...done")

elif ( "betaPlus" in process):
	print("betaPlus case not yet treated")
else:
	print("process not known. Check config file")

make_potential_from_dhfs("initial")
make_potential_from_dhfs("final")