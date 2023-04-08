import os
import numpy as np

processDict = {
  "betaMinus"  : {"DeltaZNuc" : 1,  "DeltaZAtom" : 0},
  "2betaMinus" : {"DeltaZNuc" : 2,  "DeltaZAtom" : 0},
  "betaPlus"   : {"DeltaZNuc" : -1,  "DeltaZAtom" : 0},
  "2betaPlus"  : {"DeltaZNuc" : -2,  "DeltaZAtom" : 0},
  "EC"         : {"DeltaZNuc" : -1,  "DeltaZAtom" : 1},
  "2EC"        : {"DeltaZNuc" : -2,  "DeltaZAtom" : 2},
  "ECBetaPlus" : {"DeltaZNuc" : -2,  "DeltaZAtom" : 1}
};

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
		potential_no_exchange = rv_nuc+rv_el
		
		np.savetxt(f'potential_with_Latter_{nuc}.dat', np.c_[radius, potential_with_Latter], fmt='%15.8e')
		np.savetxt(f'potential_no_Latter_{nuc}.dat', np.c_[radius, potential_no_Latter], fmt='%15.8e')
		np.savetxt(f'potential_no_exchange_{nuc}.dat', np.c_[radius, potential_no_exchange], fmt='%15.8e')


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
				line = "potentialType= "+config["DIRAC"]["potentialType"]+"\n"
			if "potentialRepository= " in line:
				line = "potentialRepository= ./"+"\n"
			if "potentialFileName= " in line:
				line = f'potentialFileName= {str(config["DIRAC"]["potentialFileName"])} \n'

			if "wfFileNameSeed= " in line:
				line = "wfFileNameSeed= " + config["DIRAC"]["wfFileNameSeed"] + "\n"
			if "outputDirectory= " in line:
				line = "outputDirectory= " + config["DIRAC"]["outputDirectory"]

			if "applyScreening= " in line:
				line = "applyScreening= "+ config["DIRAC"]["applyScreening"]+"\n"
			
			if "minimumEnergy= " in line:
				line = "minimumEnergy= "+ config["DIRAC"]["minimumEnergy"]+"\n"
			if "maximumEnergy= " in line:
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

