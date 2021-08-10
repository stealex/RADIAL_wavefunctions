#include "config_utility.hh"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <vector>
#include <math.h>

config_utility *config_utility::fInstance = 0;

config_utility *config_utility::GetInstance(){
    if (!fInstance) fInstance = new config_utility();
    return fInstance;
}

config_utility::config_utility():
processName(""),
nucleusName(""),
wfType(""),
potentialType(""),
potentialRepository(""),
potentialFileName(""),
zParent(0.),
aParent(0.),
maximumRadius(0.),
minimumRadius(0.),
minimumEnergy(0.),
maximumEnergy(0.),
nEnergyPoints(0),
nRadialPoints(0),
maxPrincipalQN(0),
minPrincipalQN(0),
applyScreening(false),
writeWF(false),
nuclearRadius(0.)
{
  kBounds[0]=0;
  kBounds[1]=0;
}

config_utility::~config_utility(){

}

void config_utility::Initialize(const std::string &configFileName){
    std::ifstream configFile(configFileName.data());
    if (!configFile.is_open()){
        std::cerr << "Could not open config file: " << configFileName.data() << "\n";
    }
    std::cout << "Reading from configuration file..." << std::endl;

    std::string line;
    while(getline(configFile, line)){

        // if line starts with "#", it's a comment. Ignore
        if (line.empty() || line[0] == '#')
            continue;
        std::vector<std::string> tokens;

        std::stringstream tmp(line);
        std::string interm;

        // Tokenize w.r.t. ' '
        while (std::getline(tmp, interm, ' '))
            tokens.push_back(interm);

        // now really read the configuration
        if (tokens[0].find("processName") != std::string::npos){
            processName = tokens[1];
            tokens.clear();

            std::cout << "Process Name = " << fInstance->processName.data() << std::endl;
            continue;
        }
        else if (tokens[0].find("nucleusName") != std::string::npos){
            nucleusName = tokens[1];
            tokens.clear();
            std::cout << "Nucleus Name = " << fInstance->nucleusName.data() << std::endl;
            continue;
        }
        else if (tokens[0].find("zParent") != std::string::npos){
            zParent = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("aParent") != std::string::npos){
            aParent = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("radiusBounds") != std::string::npos){
            minimumRadius = std::atof(tokens[1].data());
            maximumRadius = std::atof(tokens[2].data());
            tokens.clear();
            std::cout << "Radius bounds for computations: min " << fInstance->minimumRadius << " fm" <<
            " max " << fInstance->maximumRadius << " fm" <<std::endl;
            continue;
        }
        else if (tokens[0].find("nRadialPoints") != std::string::npos){
            nRadialPoints = std::atoi(tokens[1].data());
            tokens.clear();
            std::cout << "Number of radial points = " << fInstance->nRadialPoints << std::endl;
            continue;
        }
        else if (tokens[0].find("wavefunctionsType") != std::string::npos){
            wfType = tokens[1];
            tokens.clear();
            std::cout << "Will compute " << fInstance->wfType.data() << " wavefunctions" <<std::endl;
            continue;
        }
        else if (tokens[0].find("potentialType") != std::string::npos){
            potentialType = tokens[1];
            tokens.clear();
            std::cout << "Will use " << fInstance->potentialType.data() << " potential" <<std::endl;
            continue;
        }
        else if (tokens[0].find("potentialRepository") != std::string::npos){
          potentialRepository = tokens[1];
          tokens.clear();
          std::cout << "Potential repository: " << fInstance->potentialRepository.data() << std::endl;
          continue;
        }
        else if (tokens[0].find("potentialFileName") != std::string::npos){
          potentialFileName = tokens[1];
          tokens.clear();
          std::cout << "Potential prefix: " << fInstance->potentialFileName.data() << std::endl;
          continue;
        }
        else if (tokens[0].find("applyScreening") != std::string::npos){
            applyScreening = std::atoi(tokens[1].data());
            tokens.clear();
            std::cout << "Will ";
            if (applyScreening)
                std::cout << "";
            else std::cout << "NOT ";
            std::cout << "apply screening correction" <<std::endl;
            continue;
        }
        else if (tokens[0].find("minimumEnergy") != std::string::npos){
            minimumEnergy = std::atof(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("maximumEnergy") != std::string::npos){
            maximumEnergy = std::atof(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("nEnergyPoints") != std::string::npos){
            nEnergyPoints = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("writeWF") != std::string::npos){
            writeWF = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("kBounds") != std::string::npos){
            kBounds[0] = std::atoi(tokens[1].data());
            kBounds[1] = std::atoi(tokens[2].data());
            tokens.clear();
            std::cout << "Will compute for k between " << fInstance->kBounds[0] <<
             " and " << fInstance->kBounds[1] << std::endl;
            continue;
        }
        else if (tokens[0].find("maxPrincipalQN") != std::string::npos){
            maxPrincipalQN = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else if (tokens[0].find("minPrincipalQN") != std::string::npos){
            minPrincipalQN = std::atoi(tokens[1].data());
            tokens.clear();
            continue;
        }
        else{
            std::cout << "Unknown parameter: "<< tokens[0] << ". Aborting..." << std::endl;
            tokens.clear();
            exit(0);
        }
    }

    std::cout << "Parent nucleus: Z = " << zParent <<
        " A = " << aParent << std::endl;
    if (processName.find("beta") != std::string::npos){
        if (wfType.find("bound") != std::string::npos){
            std::cerr << "ERROR! Wrong configuration. Bound states only for EC" << std::endl;
            exit(1);
        }
        std::cout << "Minimum energy = " << minimumEnergy << "; Maximal energy = " << maximumEnergy << " MeV" << std::endl;
        std::cout << "Number of energy points = " << nEnergyPoints << std::endl;
    }
    if (writeWF){
        if (wfType.find("beta") != std::string::npos){
            std::cout << "Will write radial wave-functions at each energy point" << std::endl;
        }
    }
    if (wfType.find("bound") != std::string::npos){
        std::cout << "Will compute for n from " << minPrincipalQN << " up to " << maxPrincipalQN << std::endl;
    }
    configFile.close();
    std::cout << "... done reading configuration" << std::endl << std::endl;

    nuclearRadius = 1.2 * std::pow(aParent, 1./3.);
}