import os
import numpy as np


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
