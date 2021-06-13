#include "radial_interface.hh"
#include "math_tools.hh"
#include "physical_constants.hh"
#include "config_utility.hh"
#include "spline.h"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <vector>
#include <math.h>


#define splerr splerr_
#define vint vint_
#define delinf delinf_
#define sgrid sgrid_
#define radwf radwf_
#define dfree dfree_
#define dbound dbound_

#define NDIM 25000

extern "C" {
    void splerr(double *x, double *y, double *s1, double *sn, double *err, int *n, int *IWR);
    void vint(double *x, double *y, int *n);
    void delinf(double *hedel);
    void sgrid(double *r, double *dr, double *rn, double *r2, double *drn, int *n, int *nmax, int *ier);

    void dfree(double *e, double *eps, double *phase, int *k, int *irwf);
    void dbound(double *e, double *eps, int *n, int *k);

    extern struct{
        double RAD[NDIM];
        double P[NDIM];
        double Q[NDIM];
        int NGP;
        int ILAST;
        int IER;
    }radwf_;
}

using namespace physical_constants;

radial_interface *radial_interface::fInstance = 0;

radial_interface *radial_interface::GetInstance(){
    if (!fInstance) fInstance = new radial_interface();
    return fInstance;
}

radial_interface::radial_interface(){
    fConfig = config_utility::GetInstance();
}

radial_interface::~radial_interface(){
}


void radial_interface::Initialize(){
    std::cout << "Initialization phase " << std::endl;

    math_tools *mathInstance = math_tools::GetInstance();
    if (fConfig->wfType.find("scattering") != std::string::npos){
        mathInstance->ComputeEnergyPoints(
            fConfig->nEnergyPoints, fConfig->maximalEnergy, energyPoints);
    }

    // potential
    double zPotential = fConfig->zParent;
    if (fConfig->processName.find("betaMinus") != std::string::npos){
        zPotential++;
        if (fConfig->processName.find("2") != std::string::npos)
            zPotential++;
        zPotential = -1*zPotential;
    }
    else if (fConfig->processName.find("betaPlus") != std::string::npos){
        zPotential--;
        if (fConfig->processName.find("2") != std::string::npos)
            zPotential--;
    }
    else zPotential = -1. * zPotential; //EC, so electron is bound in nucleus field.

    ObtainPotential(zPotential);

    //screening
    if (fConfig->applyScreening) ApplyScreening();

    // pass potential to fortran
    double s1 = 0., sn=0.;
    double err = 0.;
    int iwr = 0;
    splerr(rValues.data(), rvValues.data(), &s1, &sn, &err, &fConfig->nRadialPoints, &iwr);
    std::cout << "Maximal error from interpolation: " << err << " (relative)" << std::endl;

    vint(rValues.data(), rvValues.data(), &fConfig->nRadialPoints);

    // High-energy limit of the dirac inner phase shift
    double hedel = 0.0;
    delinf(&hedel);
    std::cout << "Delta_infinity = " << hedel << std::endl;
}

void radial_interface::ObtainPotential(double zPotential){

    if (fConfig->potentialType.find("ChargedSphere") == std::string::npos){
        std::cout << "ERROR! only ChargedSphere potential supported currently" << std::endl;
        exit(1);
    }

    double rMinLog = std::log10(fConfig->minimumRadius);
    double rMaxLog = std::log10(fConfig->maximumRadius);
    double dr = (rMaxLog - rMinLog)/(fConfig->nRadialPoints -1.);

    math_tools *mathInstance = math_tools::GetInstance();
    rValues.clear();
    potValues.clear();
    rvValues.clear();

    rValues.push_back(0.);
    potValues.push_back(0.);
    rvValues.push_back(0.);

    for (int i = 1; i < fConfig->nRadialPoints; i++){
        rValues.push_back(std::pow(10., rMinLog + i*dr)/a0);
        potValues.push_back(mathInstance->ChargedSpherePot(zPotential, fConfig->nuclearRadius, a0*rValues[i])/e0);
        rvValues.push_back(rValues[i]*potValues[i]);
    }
}

void radial_interface::ApplyScreening(int nPoints){
    math_tools *mathInstance = math_tools::GetInstance();

    double bScreen = 0.8853 * a0 * std::pow(fConfig->zParent, -1./3.);
    double xMaxScreening = fConfig->maximumRadius / bScreen;
    mathInstance->ComputeScreeningFunction(nPoints, xMaxScreening);

    std::vector<double> screenRvalues = mathInstance->screenXvals;
    std::vector<double> phiValues = mathInstance->screenPhiVals;

    for (size_t i = 0; i < screenRvalues.size(); i++){
        screenRvalues[i] = screenRvalues[i]*bScreen/a0; // atomic units
    }

    // spline interpolation of the screening
    tk::spline s(screenRvalues, phiValues);

    for (size_t i = 0; i < rValues.size(); i++){
        double phi = s(rValues[i]);

        // correction on rV based on process
        if (fConfig->processName.find("2beta") != std::string::npos){
            // double beta... so the correction is applied in the same way
            rvValues[i] = (rvValues[i] + 2.)*phi -2.;
        }
        else if (fConfig->processName.find("beta") != std::string::npos){
            // single beta... so the correction is applied in the same way
            rvValues[i] = (rvValues[i] + 1.)*phi -1.;
        }
        else if (fConfig->processName.find("EC") != std::string::npos){
            rvValues[i] = rvValues[i] * phi;
        }
    }
}

void radial_interface::SolveDirac(){
    if (fConfig->wfType.find("scattering") != std::string::npos)
        SolveScatteringStates();
    else
        SolveBoundStates();
}

void radial_interface::SolveScatteringStates(){
    double DR0[NDIM];
    radwf.NGP = fConfig->nRadialPoints;

    // Prepare file for writting wave-functions on surface
    std::string surfaceWFFileName =
        "../Nuclei/"+fConfig->nucleusName+"/"+fConfig->nucleusName+
        "_surface_"+fConfig->wfType+"_"+fConfig->potentialType + "_screening" +
        std::to_string(fConfig->applyScreening) +
        "_" + fConfig->processName + ".dat";
    std::ofstream surfaceWFFile;

    surfaceWFFile.open(surfaceWFFileName.data(), std::ofstream::trunc);
    if (!surfaceWFFile.is_open()){
        std::cout << "ERROR! Could not open surface wave-function file. " <<
            "Check if directory exists: " << surfaceWFFileName.data() << std::endl;
        exit(1);
    }

    surfaceWFFile.width(7);
    surfaceWFFile << "E ";
    surfaceWFFile.width(7);
    surfaceWFFile << " ";

    std::vector<std::string> components = {"g","f"};
    for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++){
        if (iK == 0) continue;
        for (int iComp = 0; iComp < 2; iComp++){
            surfaceWFFile.width(10);
            surfaceWFFile<< "Re(" << components[iComp].data() << (iK < 0 ? "m" : "p" )<< abs(iK)<< ")";
            surfaceWFFile.width(10);
            surfaceWFFile<< "Im(" << components[iComp]<< (iK < 0 ? "m" : "p") << abs(iK) << ")";
        }
    }
    surfaceWFFile << std::endl;

    std::vector<double> surfWF;

    for (size_t iEnergy = 0; iEnergy < energyPoints.size(); iEnergy++){
        double e = energyPoints[iEnergy]/e0;
        if ((iEnergy+1)%10 == 0 || iEnergy == 0){
            std::cout << "  computing for energy point " << iEnergy << "; e = " << e*e0 << " MeV..." << std::endl;
        }
        // setting up the run
        double momentum = std::sqrt(energyPoints[iEnergy]*(energyPoints[iEnergy] + 2.*electronMass))/e0;
        double norm = std::sqrt(energyPoints[iEnergy]/(2.*(energyPoints[iEnergy] + electronMass)));

        double waveLength = 2.0*pi/std::sqrt(e*(2.0 + e*fineStructure*fineStructure));
        double drn = waveLength/40.0;
        double rn = drn*(radwf_.NGP - 300.);
        double r2 = 1.E-7;

        int iERR;
        int nMaxPoints = NDIM;
        // set-up grid
        sgrid(radwf.RAD, DR0, &rn, &r2, &drn, &radwf.NGP, &nMaxPoints, &iERR);
        if (iERR > 0){
            std::cout << "Problem in SGRID. iERR = " << iERR << std::endl;
        }

        // looping over requested K values
        for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++){
            if (iK == 0) continue;
//            std::cout << "   ... k = " << iK << std::endl;
            int kValue = iK;
            double eps = 1.E-15;
            double phase;
            int irwf = 1;

            dfree(&e, &eps, &phase, &kValue, &irwf);

            std::vector<double> pqSurf(0);
            FindPQSurface(pqSurf);

            double factorRe = normalUnits * std::cos(phase) * norm/
                (fConfig->nuclearRadius/a0 * momentum);
            double factorIm = normalUnits * std::sin(phase) * norm/
                (fConfig->nuclearRadius/a0 * momentum);

            surfWF.push_back(factorRe*pqSurf[0]); //Re(g)
            surfWF.push_back(factorIm*pqSurf[0]); //Im(g)
            surfWF.push_back(factorRe*pqSurf[1]); //Re(f)
            surfWF.push_back(factorIm*pqSurf[1]); //Im(f)

            if (fConfig->writeWF == 0) continue;

            std::string radFileName = "../Nuclei/" + fConfig->nucleusName + "/"+
                fConfig->nucleusName + "_" + fConfig->wfType + "_radial_" + fConfig->potentialType +
                "_screening"+ std::to_string(fConfig->applyScreening) +
                "_" + fConfig->processName+
                "_iE" + std::to_string(iEnergy) + "_k" + std::to_string(iK) + ".dat";
            std::ofstream radFile;
            radFile.open(radFileName.data(), std::ofstream::trunc);

            radFile << std::scientific << e*e0 << std::endl;
            radFile << std::scientific << phase << std::endl;

            int irMax = -1;
            for (int i = radwf.NGP-1; i >= 0; i--){
                if (abs(radwf.P[i])  > 1.E-35){
                    irMax = i;
                    break;
                }
            }

            WriteBoundWFFile(irMax, radFile);
            radFile.close();
        }

        WriteSurfWFLine(energyPoints[iEnergy] + electronMass, surfWF, surfaceWFFile);
        surfWF.clear();
    }

    surfaceWFFile.close();
}
void radial_interface::WriteSurfWFLine(double e, std::vector<double> vals, std::ofstream &file){
    file.precision(6);
    file.width(14);
    file << e;

    for (size_t i = 0; i < vals.size(); i++){
        file.precision(6);
        file.width(14);
        file << std::scientific;
        file <<  vals[i];
    }
    file << std::endl;
}

void radial_interface::FindPQSurface(std::vector<double> &pq){
// find surface point
    int rSurfBounds[2] = {-1, -1};

    for (int iR = 1; iR < NDIM; iR++){
        if (radwf.RAD[iR]*a0 >= fConfig->nuclearRadius){
            rSurfBounds[0] = iR -1;
            break;
        }
    }

    for (int iR = 1; iR < NDIM; iR++){
        if (radwf.RAD[iR]*a0 > fConfig->nuclearRadius){
            rSurfBounds[1] = iR ;
            break;
        }
    }

    if (rSurfBounds[0] < 0 || rSurfBounds[1] < 0){
        std::cout << "Problem finding surface" << std::endl;
    }

    double slopeP = (radwf.P[rSurfBounds[1]]-radwf.P[rSurfBounds[0]])/(radwf.RAD[rSurfBounds[1]]-radwf.RAD[rSurfBounds[0]]);
    double slopeQ = (radwf.Q[rSurfBounds[1]]-radwf.Q[rSurfBounds[0]])/(radwf.RAD[rSurfBounds[1]]-radwf.RAD[rSurfBounds[0]]);

    double dr = fConfig->nuclearRadius/a0 - radwf.RAD[rSurfBounds[0]];
    double pSurf = slopeP * dr + radwf.P[rSurfBounds[0]];
    double qSurf = slopeQ * dr + radwf.Q[rSurfBounds[0]];

    pq.push_back(pSurf);
    pq.push_back(qSurf);
}

void radial_interface::SolveBoundStates(){
    //setup for bound states computation
    radwf.NGP = fConfig->nRadialPoints;
    double rn = 10.0; //atomic units
    double DR0[NDIM];
    double r2 = 1.0E-6; //atomic units. Must be larger than 1E-8
    double drn = 1.0;
    int iERR;
    int nMaxPoints = NDIM;
    double eps = 1.E-12;
    // grid
    sgrid(radwf.RAD, DR0, &rn, &r2, &drn, &radwf.NGP, &nMaxPoints, &iERR);
    if (iERR > 0){
        std::cout << "Problem with SGRID. iERR = " << iERR << std::endl;
    }

    //actual computation
    for (int iN = 1; iN <= fConfig->maxPrincipalQN; iN++){
        std::cout << "Computing for n = " << iN <<"\n";
        int nValue = iN;


        for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++){
            if (iK == 0 || iK < -iN || iK > iN -1) continue;
            int kValue = iK;
            std::string boundWFFileName =
                "../Nuclei/" + fConfig->nucleusName + "/" + fConfig->nucleusName +
                "_" + fConfig->wfType + "_" + fConfig->potentialType + "_screening" +
                std::to_string(fConfig->applyScreening) + "_" + fConfig->processName +
                "_n"+std::to_string(iN) + "_k" + std::to_string(iK) + ".dat";

            std::ofstream boundWFFile;
            boundWFFile.open(boundWFFileName, std::ofstream::trunc);

            std::cout << " ... k = " << iK << "\n";

            double e = -1*std::pow(fConfig->zParent, 2.)/(2.0 * iN * iN);
            std::cout << "   Trial energy " << std::scientific << e*e0 << "MeV\n";
            dbound(&e, &eps, &nValue, &kValue);
            std::cout << "   Bound energy " << std::scientific << e*e0 << "MeV\n" ;

            boundWFFile << std::scientific << e*e0 << std::endl;

            int iMaxR = -1;
            for (int i = radwf.NGP-1; i >= 0; i--){
                if (abs(radwf.P[i]) > 1.E-35){
                    iMaxR = i;
                    break;
                }
            }

            WriteBoundWFFile(iMaxR, boundWFFile);
            boundWFFile.close();
        }
    }
}

void radial_interface::WriteBoundWFFile(int irMax, std::ofstream &file){

    file.width(7);
    file << "R ";
    file.width(7);
    file << " ";

    std::vector<std::string> components = {"P","Q"};
    for (int iComp = 0; iComp < 2; iComp++){
        file.width(10);
        file << components[iComp].data();
    }
    file << std::endl;

    // write the values of R P and Q
    for (int iR = 0; iR < irMax; iR++){
        file.precision(6);
        file.width(14);
        file << radwf.RAD[iR];

        file.precision(6);
        file.width(14);
        file << std::scientific;
        file <<  radwf.P[iR];
        file.precision(6);
        file.width(14);
        file << std::scientific;
        file << radwf.Q[iR];

        file << std::endl;
    }
}
