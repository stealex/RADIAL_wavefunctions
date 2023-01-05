import configparser
import subprocess
import numpy as np
import os
import shutil
from mendeleev import element
import re
from python_util import dhfs_util

#
#  Helper functions
#
#
# Script main
#
def main():
  project_base_dir = os.getcwd()
  print("base dir = ", project_base_dir)

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
  dhfs_util.move_dhfs_output(z_initial, "initial")

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

    dhfs_util.move_dhfs_output(z_final, "final")
    print("...done")

  elif ( "betaPlus" in process):
    print("betaPlus case not yet treated")
  else:
    print("process not known. Check config file")

  dhfs_util.make_potential_from_dhfs("initial")
  dhfs_util.make_potential_from_dhfs("final")

  # create config file for RADIAL initial nucleus
  initial_potential_file = ""
  initial_wfFileName_seed = ""
  final_potential_file = ""
  final_wfFileName_seed = ""
  if config["DHFS"]["LATTER_TAIL"] == "1":
    initial_potential_file = "potential_with_Latter_initial.dat"
    initial_wfFileName_seed = "initial_DHFS_LATTER_TAIL_1"
    final_potential_file = "potential_with_Latter_final.dat"
    final_wfFileName_seed = "final_DHFS_LATTER_TAIL_1"
  elif config["DHFS"]["LATTER_TAIL"] == "0":
    initial_potential_file = "potential_no_Latter_initial.dat"
    initial_wfFileName_seed = "initial_DHFS_LATTER_TAIL_0"
    final_potential_file = "potential_with_Latter_final.dat"
    final_wfFileName_seed = "final_DHFS_LATTER_TAIL_0"
  else:
    print("WRONG configuration. LATTER_TAIL can be 0 or 1")
    exit(1)

  config["DIRAC"]["outputDirectory"] = "."
  config["DIRAC"]["wavefunctionsType"] = "bound"
  config["DIRAC"]["processName"] = "EC"
  config["DIRAC"]["potentialFileName"] = initial_potential_file
  config["DIRAC"]["wfFileNameSeed"] = initial_wfFileName_seed

  dhfs_util.create_dirac_config("wavefunctions_initial_with_Latter.conf", config)
  subprocess.run(["./Radial_WF", "wavefunctions_initial_with_Latter.conf"])

  # create config file for RADIAL final nucleus

  ## bound states
  config["DIRAC"]["potentialFileName"] = final_potential_file
  config["DIRAC"]["wfFileNameSeed"] = final_wfFileName_seed

  dhfs_util.create_dirac_config("wavefunctions_final_with_Latter.conf", config)
  subprocess.run(["./Radial_WF", "wavefunctions_final_with_Latter.conf"])

  ## scattering states
  config["DIRAC"]["wavefunctionsType"] = "scattering"
  config["DIRAC"]["processName"] = "betaMinus"
  config["DIRAC"]["potentialFileName"] = final_potential_file
  config["DIRAC"]["wfFileNameSeed"] = final_wfFileName_seed

  dhfs_util.create_dirac_config("wavefunctions_final_with_Latter_scattering.conf", config)
  subprocess.run(["./Radial_WF", "wavefunctions_final_with_Latter_scattering.conf"])


if __name__ == "__main__":
  main()