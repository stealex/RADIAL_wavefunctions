#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <string>

#include "radial_interface.hh"
#include "config_utility.hh"

void Usage(){
    printf("./bin/main <config_file>\n");
}

int main(int argc, char **argv){
    if (argc != 2){
        Usage();
    }

    // initialize configuration
    config_utility *config = config_utility::GetInstance();
    std::string configFileName = argv[1];
    config->Initialize(configFileName);

    // initialize radial interface
    radial_interface *interface = radial_interface::GetInstance();
    interface->Initialize();

    // solve Dirac equation. Magic happens in radial_interface.cc
    interface->SolveDirac();
}