import configparser
import subprocess
import numpy as np
import os
import shutil
from mendeleev import element
import re


#
#  Helper functions
#
def move_dhfs_output(z, nuc):
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
		
		np.savetxt(f'potential_with_Latter_{nuc}.dat', np.c_[radius, potential_with_Latter], fmt='%15.8e')
		np.savetxt(f'potential_no_Latter_{nuc}.dat', np.c_[radius, potential_no_Latter], fmt='%15.8e')

def create_dirac_config(out_file_name, config):
	with open("wavefunctions_initial_tmp.conf") as config_file_tmp:
		lines = config_file_tmp.readlines()
		config_file = open(out_file_name, 'w')
		for line in lines:
			if "processName= " in line:
				line = "processName= " + config["DIRAC"]["processName"]+"\n"
			
			if "nucleusName= " in line:
				line = "nucleusName= " + config["DIRAC"]["nucleusName"]+"\n"
			if "zParent= " in line:
				line = "zParent= " + config["DIRAC"]["zInitial"]+"\n"
			if "aParent= " in line:
				line = "aParent= " + config["DIRAC"]["aInitial"]+"\n"
			
			if "radiusBounds= " in line:
				line = "radiusBounds= " + config["DIRAC"]["radiusBounds"]+"\n"
			if "nRadialPoints= " in line:
				line = "nRadialPoints= " + config["DIRAC"]["nRadialPoints"]+"\n"
			
			if "wavefunctionsType= " in line:
				line = "wavefunctionsType= "+ config["DIRAC"]["wavefunctionsType"]+"\n"

			if "potentialType= " in line:
				line = "potentialType= FromFile"+"\n"
			if "potentialRepository= " in line:
				line = "potentialRepository= ./"+"\n"
			if "potentialFileName= " in line:
				line = f'potentialFileName= {str(config["DIRAC"]["potentialFileName"])} \n'

			if "wfFileNameSeed= " in line:
				line = "wfFileNameSeed= " + config["DIRAC"]["wfFileNameSeed"] + "\n"
			if "outputDirectory= " in line:
				line = "outputDirectory= " + config["DIRAC"]["outputDirectory"]

			if "applyScreening= " in line:
				line = "applyScreening= 0"+"\n"
			
			if "minimumEnergy= " in line:
				line = "minimumEnergy= "+ config["DIRAC"]["minimumEnergy"]+"\n"
			if "maxmimumEnergy= " in line:
				line = "maximumEnergy= "+ config["DIRAC"]["maximumEnergy"]+"\n"
			if "nEnergyPoints= " in line:
				line = "nEnergyPoints= " + config["DIRAC"]["nEnergyPoints"]+"\n"
			
			if "writeWF= " in line:
				line = "writeWF= "+config["DIRAC"]["writeWF"]+"\n"
			if "kBounds= " in line:
				line = "kBounds= "+config["DIRAC"]["kBounds"]+"\n"

			if "minPrincipalQN= " in line:
				line = "minPrincipalQN= "+config["DIRAC"]["minPrincipalQN"]+"\n"
			if "maxPrincipalQN= " in line:
				line = "maxPrincipalQN= "+config["DIRAC"]["maxPrincipalQN"]+"\n"

			config_file.write(line)
		config_file.close()

#
# Script main
#

project_base_dir = os.getcwd()
print("base dir = ", project_base_dir)

config = configparser.ConfigParser()
config.read("dhfs.ini")

# initial nucleus name and the process it undergoes
initial_nucleus = config["DIRAC"]["nucleusName"]
process = config["DIRAC"]["processName"]
z_initial = int(config["DIRAC"]["zInitial"])
a_initial = int(config["DIRAC"]["aInitial"])
initial_element = element(z_initial)

# create directory to store results
final_dir = initial_nucleus+"_"+process
try:
	os.mkdir(final_dir)
except OSError as error:
	print(error)

shutil.copy2(project_base_dir+"/build/DHFS", final_dir+"/DHFS")
shutil.copy2(project_base_dir+"/build/constants.mod", final_dir+"/constants.mod")
shutil.copy2(project_base_dir+"/build/Radial_WF", final_dir+"/Radial_WF")

# copy the dhfs .in file to the final directory
if int(config["DHFS"]["USE_STANDARD"]) == 1:
	repo = project_base_dir+"/"+config["DHFS"]["REPO"]
	shutil.copyfile(repo+"/"+f'Z{z_initial:0>3d}.in', final_dir+"/dhfs_initial.in")
else:
	shutil.copyfile(str(config["INPUT_FILE"]), final_dir+"/dhfs_initial.in")

shutil.copyfile(str(project_base_dir+"/wavefunctions.conf"), final_dir+"/wavefunctions_initial_tmp.conf")

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

# create config file for RADIAL initial nucleus
# TODO: get rid of process name
config["DIRAC"]["outputDirectory"] = "."
config["DIRAC"]["wavefunctionsType"] = "bound"
config["DIRAC"]["processName"] = "EC"

if config["DHFS"]["LATTER_TAIL"] == "1":
  config["DIRAC"]["potentialFileName"] = "potential_with_Latter_initial.dat"
  config["DIRAC"]["wfFileNameSeed"] = "initial_DHFS_LATTER_TAIL_1"
elif config["DHFS"]["LATTER_TAIL"] == "0":
  config["DIRAC"]["potentialFileName"] = "potential_no_Latter_initial.dat"
  config["DIRAC"]["wfFileNameSeed"] = "initial_DHFS_LATTER_TAIL_0"
else:
  print("WRONG configuration. LATTER_TAIL can be 0 or 1")
  exit(1)

create_dirac_config("wavefunctions_initial_with_Latter.conf", config)
subprocess.run(["./Radial_WF", "wavefunctions_initial_with_Latter.conf"])


# create config file for RADIAL final nucleus

## bound states
config["DIRAC"]["outputDirectory"] = "."
if config["DHFS"]["LATTER_TAIL"] == "1":
  config["DIRAC"]["potentialFileName"] = "potential_with_Latter_final.dat"
  config["DIRAC"]["wfFileNameSeed"] = "final_DHFS_LATTER_TAIL_1"
elif config["DHFS"]["LATTER_TAIL"] == "0":
  config["DIRAC"]["potentialFileName"] = "potential_no_Latter_final.dat"
  config["DIRAC"]["wfFileNameSeed"] = "final_DHFS_LATTER_TAIL_0"
else:
  print("WRONG configuration. LATTER_TAIL can be 0 or 1")
  exit(1)

create_dirac_config("wavefunctions_final_with_Latter.conf", config)
subprocess.run(["./Radial_WF", "wavefunctions_final_with_Latter.conf"])

## scattering states
config["DIRAC"]["outputDirectory"] = "."
config["DIRAC"]["wavefunctionsType"] = "scattering"
config["DIRAC"]["processName"] = "betaMinus"
if config["DHFS"]["LATTER_TAIL"] == "1":
  config["DIRAC"]["potentialFileName"] = "potential_with_Latter_final.dat"
  config["DIRAC"]["wfFileNameSeed"] = "final_DHFS_LATTER_TAIL_1"
elif config["DHFS"]["LATTER_TAIL"] == "0":
  config["DIRAC"]["potentialFileName"] = "potential_no_Latter_final.dat"
  config["DIRAC"]["wfFileNameSeed"] = "final_DHFS_LATTER_TAIL_0"
else:
  print("WRONG configuration. LATTER_TAIL can be 0 or 1")
  exit(1)

create_dirac_config("wavefunctions_final_with_Latter_scattering.conf", config)
subprocess.run(["./Radial_WF", "wavefunctions_final_with_Latter_scattering.conf"])
