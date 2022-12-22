#include "radial_interface.hh"
#include "config_utility.hh"
#include "math_tools.hh"
#include "physical_constants.hh"
#include "spline.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <vector>

#define splerr splerr_
#define vint vint_
#define delinf delinf_
#define sgrid sgrid_
#define radwf radwf_
#define dfree dfree_
#define dbound dbound_

#define NDIM 25000

extern "C" {
void splerr(double *x, double *y, double *s1, double *sn, double *err, int *n,
            int *IWR);
void vint(double *x, double *y, int *n);
void delinf(double *hedel);
void sgrid(double *r, double *dr, double *rn, double *r2, double *drn, int *n,
           int *nmax, int *ier);

void dfree(double *e, double *eps, double *phase, int *k, int *irwf);
void dbound(double *e, double *eps, int *n, int *k);

extern struct {
  double RAD[NDIM];
  double P[NDIM];
  double Q[NDIM];
  int    NGP;
  int    ILAST;
  int    IER;
} radwf_;
}

using namespace physical_constants;

radial_interface *radial_interface::fInstance = 0;

radial_interface *radial_interface::GetInstance() {
  if (!fInstance)
    fInstance = new radial_interface();
  return fInstance;
}

radial_interface::radial_interface() {
  fConfig = config_utility::GetInstance();
}

radial_interface::~radial_interface() {}

void radial_interface::Initialize() {
  std::cout << "Initialization phase " << std::endl;

  math_tools *mathInstance = math_tools::GetInstance();
  if (fConfig->wfType.find("scattering") != std::string::npos) {
    mathInstance->ComputeEnergyPoints(fConfig->nEnergyPoints,
                                      fConfig->minimumEnergy, fConfig->maximumEnergy,
                                      energyPoints);
  }

  // potential
  double zPotential = fConfig->zParent;
  if (fConfig->processName.find("betaMinus") != std::string::npos) {
    zPotential++;
    if (fConfig->processName.find("2") != std::string::npos)
      zPotential++;
    zPotential = -1 * zPotential;
  } else if (fConfig->processName.find("betaPlus") != std::string::npos) {
    zPotential--;
    if (fConfig->processName.find("2") != std::string::npos)
      zPotential--;
  } else
    zPotential =
      -1. * fConfig->zParent; // EC, so electron is bound in nucleus field.

  ObtainPotential(zPotential);

  // screening
  if (fConfig->applyScreening)
    ApplyScreening();

  // pass potential to fortran
  double s1 = 0., sn = 0.;
  double err = 0.;
  int    iwr = 0;
  splerr(rValues.data(), rvValues.data(), &s1, &sn, &err,
         &fConfig->nRadialPoints, &iwr);
  std::cout << "Maximal error from interpolation: " << err << " (relative)"
            << std::endl;

  vint(rValues.data(), rvValues.data(), &fConfig->nRadialPoints);

  // High-energy limit of the dirac inner phase shift
  double hedel = 0.0;
  delinf(&hedel);
  std::cout << "Delta_infinity = " << hedel << std::endl;
}

void radial_interface::ObtainPotential(double zPotential) {

  if (fConfig->potentialType.find("ChargedSphere") == std::string::npos &&
      fConfig->potentialType.find("PointCharge") == std::string::npos && 
      fConfig->potentialType.find("FromFile") == std::string::npos) {
    std::cout
      << "ERROR! only ChargedSphere, PointCharge, FromFile potential supported currently"
      << std::endl;
    exit(1);
  }

  if (fConfig->potentialType.find("ChargedSphere") != std::string::npos ||
      fConfig->potentialType.find("PointCharge") != std::string::npos){
    
    double rNuc = 0.0;
    if (fConfig->potentialType.find("ChargedSphere") != std::string::npos) {
      rNuc = fConfig->nuclearRadius;
    }
    ComputeChargedSpherePotential(rNuc, zPotential);
  }
  else if (fConfig->potentialType.find("FromFile") != std::string::npos){
    std::string potFileName = fConfig->potentialRepository + "/" + fConfig->potentialFileName;
    // for debug.
    std::ifstream potFile(potFileName.data());
    if (!potFile.is_open()){
      std::cout
        << "ERROR! could not open potential file: " << potFileName.data()
        << std::endl;
      exit(1);
    }
    std::cout << "Reading potential and radii from file: " << potFileName.data() << std::endl;
    ReadPotentialFile(potFile); 
    potFile.close();
  }
}

void radial_interface::ComputeChargedSpherePotential(double rNuc, double zPotential){
  // Computes the potential of a charged sphere of radius rNuc. 
  // If rNuc is 0, it is the potential of a point charge
  double rMinLog = std::log10(fConfig->minimumRadius);
  double rMaxLog = std::log10(fConfig->maximumRadius);
  double dr      = (rMaxLog - rMinLog) / (fConfig->nRadialPoints - 1.);

  math_tools *mathInstance = math_tools::GetInstance();
  rValues.clear();
  potValues.clear();
  rvValues.clear();

  rValues.push_back(0.);
  potValues.push_back(0.);
  rvValues.push_back(0.);

  for (int i = 1; i < fConfig->nRadialPoints; i++) {
    rValues.push_back(std::pow(10., rMinLog + i * dr) / a0);
    potValues.push_back(
      mathInstance->ChargedSpherePot(zPotential, rNuc, a0 * rValues[i]) / e0);
    rvValues.push_back(rValues[i] * potValues[i]);
  }
}

void radial_interface::ReadPotentialFile(std::ifstream &file){
  double r, rv;
  std::vector<double> rVecTmp;
  std::vector<double> vVecTmp;
  while(file >> r >> rv){
    rVecTmp.push_back(r);
    vVecTmp.push_back(rv);
  }

  std::cout << "Interpolating potential from file" << std::endl;
  // spline interpolation of the screening
  tk::spline vInterpolated(rVecTmp, vVecTmp);

  std::cout << "Constructing radial grid and potential values" << std::endl;
  
  fConfig->minimumRadius = std::max(fConfig->minimumRadius, a0*rVecTmp[0]);
  fConfig->maximumRadius = std::min(fConfig->maximumRadius, a0*rVecTmp[rVecTmp.size() - 1]);

  double rMinLog = std::log10(fConfig->minimumRadius);
  double rMaxLog = std::log10(fConfig->maximumRadius);
  double dr      = (rMaxLog - rMinLog) / (fConfig->nRadialPoints - 1.);

  rValues.push_back(0.);
  for (int i = 1; i < fConfig->nRadialPoints; i++) {
    rValues.push_back(std::pow(10., rMinLog + i * dr) / a0);
  }

  rvValues.push_back(0.);
  for (size_t i = 1; i < rValues.size(); i++) {
    rvValues.push_back(vInterpolated(rValues[i]));
  }

  std::cout << "Done reading potential from file" << std::endl;
}

void radial_interface::ApplyScreening(int nPoints) {
  math_tools *mathInstance = math_tools::GetInstance();

  double bScreen       = 0.8853 * a0 * std::pow(fConfig->zParent, -1. / 3.);
  double xMaxScreening = fConfig->maximumRadius / bScreen;
  mathInstance->ComputeScreeningFunction(nPoints, xMaxScreening);

  std::vector<double> screenRvalues = mathInstance->screenXvals;
  std::vector<double> phiValues     = mathInstance->screenPhiVals;

  for (size_t i = 0; i < screenRvalues.size(); i++) {
    screenRvalues[i] = screenRvalues[i] * bScreen / a0; // atomic units
  }

  // spline interpolation of the screening
  tk::spline s(screenRvalues, phiValues);

  for (size_t i = 0; i < rValues.size(); i++) {
    double phi = s(rValues[i]);

    // correction on rV based on process
    if (fConfig->processName.find("2beta") != std::string::npos) {
      // double beta... so the correction is applied in the same way
      rvValues[i] = (rvValues[i] + 2.) * phi - 2.;
    } else if (fConfig->processName.find("beta") != std::string::npos) {
      // single beta... so the correction is applied in the same way
      rvValues[i] = (rvValues[i] + 1.) * phi - 1.;
    } else if (fConfig->processName.find("EC") != std::string::npos) {
      rvValues[i] = rvValues[i] * phi;
    }
  }
}

void radial_interface::SolveDirac() {
  if (fConfig->wfType.find("scattering") != std::string::npos)
    SolveScatteringStates();
  else
    SolveBoundStates();
}

std::string radial_interface::CreateFileName(std::string placeOfComputation, int iEnergy = 0, int iN = 0, int iK = 0){
  std::string nameSeed = fConfig->wfFileNameSeed == "" ? fConfig->potentialType : fConfig->wfFileNameSeed; 
  std::string outFileName = "";
  std::string outputDirectory = fConfig->outputDirectory == ""? "../Nuclei" : fConfig->outputDirectory + "/";
  outFileName += outputDirectory;

  if (fConfig->outputDirectory == "")
  {
    outFileName += fConfig->nucleusName +
      "_" + placeOfComputation + "_" + fConfig->wfType + "_" + nameSeed +
      "_screening" + std::to_string(fConfig->applyScreening) + "_" +
      fConfig->processName;
  }
  else
  {
    outFileName += fConfig->wfType +
      "_" + placeOfComputation + "_" + nameSeed;
  }

  if (placeOfComputation == "radial")
  {
    if (fConfig->wfType == "scattering")
    {
      outFileName += "_iE" + std::to_string(iEnergy) + "_k" + std::to_string(iK) + ".dat";
    }
    else
    {
      outFileName += "_n" + std::to_string(iN) + "_k" + std::to_string(iK) + ".dat";
    }
  }
  else
  {
    outFileName += ".dat";
  }
  return outFileName;
}
void radial_interface::SolveScatteringStates() {
  double DR0[NDIM];
  radwf.NGP = fConfig->nRadialPoints;

  // Prepare file for writting wave-functions on surface
  std::string nameSeed = fConfig->wfFileNameSeed == "" ? fConfig->potentialType : fConfig->wfFileNameSeed;
  std::string surfaceWFFileName = CreateFileName("surface");
  
  // if (fConfig->outputDirectory == ""){
  //   surfaceWFFileName =
  //     "../Nuclei/" + fConfig->nucleusName + "/" + fConfig->nucleusName +
  //     "_surface_" + fConfig->wfType + "_" + nameSeed +
  //     "_screening" + std::to_string(fConfig->applyScreening) + "_" +
  //     fConfig->processName + ".dat";
  // }
  // else{
  //   surfaceWFFileName = fConfig->outputDirectory + "/" + fConfig->wfType +
  //     "_surface_" + nameSeed + ".dat";
  // }
  std::ofstream surfaceWFFile;

  surfaceWFFile.open(surfaceWFFileName.data(), std::ofstream::trunc);
  if (!surfaceWFFile.is_open()) {
    std::cout << "ERROR! Could not open surface wave-function file. "
              << "Check if directory exists: " << surfaceWFFileName.data()
              << std::endl;
    exit(1);
  }

  surfaceWFFile.width(7);
  surfaceWFFile << "E ";
  surfaceWFFile.width(7);
  surfaceWFFile << " ";

  std::vector<std::string> components = {"g", "f"};
  for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++) {
    if (iK == 0)
      continue;
    for (int iComp = 0; iComp < 2; iComp++) {
      surfaceWFFile.width(10);
      surfaceWFFile << "Re(" << components[iComp].data() << (iK < 0 ? "m" : "p")
                    << abs(iK) << ")";
      surfaceWFFile.width(10);
      surfaceWFFile << "Im(" << components[iComp] << (iK < 0 ? "m" : "p")
                    << abs(iK) << ")";
    }
  }
  surfaceWFFile << std::endl;

  std::vector<double> surfWF;

  for (size_t iEnergy = 0; iEnergy < energyPoints.size(); iEnergy++) {
    double e = energyPoints[iEnergy] / e0;
    if ((iEnergy + 1) % 10 == 0 || iEnergy == 0) {
      std::cout << "  computing for energy point " << iEnergy
                << "; e = " << e * e0 << " MeV..." << std::endl;
    }
    // setting up the run
    double momentum = std::sqrt(energyPoints[iEnergy] *
                                (energyPoints[iEnergy] + 2. * electronMass)) /
                      e0;
    double norm = std::sqrt((energyPoints[iEnergy] + 2. * electronMass) /
                            (2. * (energyPoints[iEnergy] + electronMass)));

    double waveLength =
      2.0 * pi / std::sqrt(e * (2.0 + e * fineStructure * fineStructure));
    // double drn = waveLength / 40.0;
    // double rn  = drn * (radwf.NGP - 300.);
    double rn = fConfig->maximumRadius/a0;
    double drn = rn/(radwf.NGP - 300);
    double r2  = 1.E-7;

    int iERR;
    int nMaxPoints = NDIM;
    // set-up grid
    sgrid(radwf.RAD, DR0, &rn, &r2, &drn, &radwf.NGP, &nMaxPoints, &iERR);
    if (iERR > 0) {
      std::cout << "Problem in SGRID. iERR = " << iERR << std::endl;
    }

    // looping over requested K values
    for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++) {
      if (iK == 0)
        continue;
      int    kValue = iK;
      double eps    = 1.E-15;
      double phase;
      int    irwf = 1;

      dfree(&e, &eps, &phase, &kValue, &irwf);
      std::vector<double> pqSurf(0);
      FindPQSurface(pqSurf);

      double factorRe = normalUnits * std::cos(phase) * norm /
                        (fConfig->nuclearRadius / a0 * momentum);
      double factorIm = normalUnits * std::sin(phase) * norm /
                        (fConfig->nuclearRadius / a0 * momentum);

      surfWF.push_back(factorRe * pqSurf[0]); // Re(g)
      surfWF.push_back(factorIm * pqSurf[0]); // Im(g)
      surfWF.push_back(factorRe * pqSurf[1]); // Re(f)
      surfWF.push_back(factorIm * pqSurf[1]); // Im(f)

      if (fConfig->writeWF == 0)
        continue;

      std::string radFileName = CreateFileName("radial", iEnergy, 0, iK);
      
      // if (fConfig->outputDirectory == ""){
      //   radFileName = "../Nuclei/" + fConfig->nucleusName + "/" + fConfig->nucleusName + "_" +
      //   fConfig->wfType + "_radial_" + nameSeed + "_screening" +
      //   std::to_string(fConfig->applyScreening) + "_" + fConfig->processName +
      //   "_iE" + std::to_string(iEnergy) + "_k" + std::to_string(iK) + ".dat";
      // }
      // else{
      //   radFileName = fConfig->outputDirectory + "/" + fConfig->nucleusName + "_" +
      //   fConfig->wfType + "_radial_" + nameSeed + "_screening" +
      //   std::to_string(fConfig->applyScreening) + "_" + fConfig->processName +
      //   "_iE" + std::to_string(iEnergy) + "_k" + std::to_string(iK) + ".dat";
      // }
      int irMax = radwf.NGP - 1;
      // for (int i = radwf.NGP - 1; i >= 0; i--) {
      //   if (abs(radwf.P[i]) > 1.E-35) {
      //     irMax = i;
      //     break;
      //   }
      // }

      WriteRadialWFFile(irMax, e*e0, radFileName, phase);
    }

    WriteSurfWFLine(energyPoints[iEnergy] + electronMass, surfWF,
                    surfaceWFFile);
    surfWF.clear();
  }

  surfaceWFFile.close();
}
void radial_interface::WriteSurfWFLine(double e, std::vector<double> vals,
                                       std::ofstream &file) {
  if (fConfig->wfType.find("scattering") != std::string::npos) {
    file.precision(6);
    file.width(14);
    file << std::scientific;
    file << e;
  }

  for (size_t i = 0; i < vals.size(); i++) {
    file.precision(6);
    file.width(14);
    file << std::scientific;
    file << vals[i];
  }
  file << std::endl;
}

void radial_interface::FindPQSurface(std::vector<double> &pq) {
  // find surface point
  int rSurfBounds[2] = {-1, -1};

  for (int iR = 1; iR < NDIM; iR++) {
    if (radwf.RAD[iR] * a0 >= fConfig->nuclearRadius) {
      rSurfBounds[0] = iR - 1;
      break;
    }
  }

  for (int iR = 1; iR < NDIM; iR++) {
    if (radwf.RAD[iR] * a0 > fConfig->nuclearRadius) {
      rSurfBounds[1] = iR;
      break;
    }
  }

  if (rSurfBounds[0] < 0 || rSurfBounds[1] < 0) {
    std::cout << "Problem finding surface" << std::endl;
  }

  double slopeP = (radwf.P[rSurfBounds[1]] - radwf.P[rSurfBounds[0]]) /
                  (radwf.RAD[rSurfBounds[1]] - radwf.RAD[rSurfBounds[0]]);
  double slopeQ = (radwf.Q[rSurfBounds[1]] - radwf.Q[rSurfBounds[0]]) /
                  (radwf.RAD[rSurfBounds[1]] - radwf.RAD[rSurfBounds[0]]);

  double dr    = fConfig->nuclearRadius / a0 - radwf.RAD[rSurfBounds[0]];
  double pSurf = slopeP * dr + radwf.P[rSurfBounds[0]];
  double qSurf = slopeQ * dr + radwf.Q[rSurfBounds[0]];

  pq.push_back(pSurf);
  pq.push_back(qSurf);
}

void radial_interface::SolveBoundStates() {
  // setup for bound states computation
  radwf.NGP = fConfig->nRadialPoints;
  double rn = fConfig->maximumRadius/a0; // atomic units
  double DR0[NDIM];
  double r2  = fConfig->minimumRadius/a0; // atomic units. Must be larger than 1E-8
  double drn = 1.0;
  int    iERR;
  int    nMaxPoints = NDIM;
  double eps        = 1.E-12;

  // prepare file for writing wave-functions on surface
  std::string nameSeed = fConfig->wfFileNameSeed == "" ? fConfig->potentialType : fConfig->wfFileNameSeed;
  std::string surfaceWFFileName = CreateFileName("surface");
  
  // if (fConfig->outputDirectory == ""){
  //   surfaceWFFileName =
  //     "../Nuclei/" + fConfig->nucleusName + "/" + fConfig->nucleusName +
  //     "_surface_" + fConfig->wfType + "_" + nameSeed +
  //     "_screening" + std::to_string(fConfig->applyScreening) + "_" +
  //     fConfig->processName + ".dat";
  // }
  // else{
  //   surfaceWFFileName = fConfig->outputDirectory + "/" + fConfig->nucleusName +
  //     "_surface_" + fConfig->wfType + "_" + nameSeed +
  //     "_screening" + std::to_string(fConfig->applyScreening) + "_" +
  //     fConfig->processName + ".dat";
  // }

  std::ofstream surfaceWFFile;
  surfaceWFFile.open(surfaceWFFileName.data(), std::ofstream::trunc);
  if (!surfaceWFFile.is_open()) {
    std::cout << "ERROR! Could not open surface wave-function file. "
              << "Check if directory exists: " << surfaceWFFileName.data()
              << std::endl;
    exit(1);
  }

  for (int iN = fConfig->minPrincipalQN; iN <= fConfig->maxPrincipalQN; iN++) {
    for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++) {
      if (iK == 0 || iK < -iN || iK > iN - 1)
        continue;
      std::string stringE = "E("+std::to_string(iN) + ","+std::to_string(iK) + ")";
      std::string stringProb = "Prob("+std::to_string(iN) + ","+std::to_string(iK)+")";

      surfaceWFFile.width(14);
      surfaceWFFile << stringE.data();
      surfaceWFFile.width(14);
      surfaceWFFile << stringProb.data();
    }
  }
  surfaceWFFile << std::endl;
  // grid
  sgrid(radwf.RAD, DR0, &rn, &r2, &drn, &radwf.NGP, &nMaxPoints, &iERR);
  if (iERR > 0) {
    std::cout << "Problem with SGRID. iERR = " << iERR << std::endl;
  }

  // actual computation
  std::vector<double> surfaceData(0);
  for (int iN = fConfig->minPrincipalQN; iN <= fConfig->maxPrincipalQN; iN++) {
    std::cout << "Computing for n = " << iN << "\n";
    int nValue = iN;

    for (int iK = fConfig->kBounds[0]; iK <= fConfig->kBounds[1]; iK++) {
      if (iK == 0 || iK < -iN || iK > iN - 1)
        continue;
      int kValue = iK;

      std::cout << " ... k = " << iK << "\n";

      // double e = -1*std::pow(fConfig->zParent, 2.)/(2.0 * iN * iN);
      double e =
        -1. * (electronMass / e0) *
        (1. -
         1. / std::sqrt(1. + std::pow(fineStructure * fConfig->zParent, 2) /
                               std::pow(iN - std::abs(iK) +
                                          std::sqrt(iK * iK -
                                                    std::pow(fineStructure *
                                                               fConfig->zParent,
                                                             2.)),
                                        2)));
      std::cout << "   Trial energy " << std::scientific << e * e0 << "MeV\n";
      dbound(&e, &eps, &nValue, &kValue);
      std::cout << "   Bound energy " << std::scientific << e * e0 << "MeV\n";

      // surface data
      std::vector<double> pqSurf(0);
      FindPQSurface(pqSurf);

      double probabilitySurface =
        std::pow(hc / (electronMass * a0), 3) *
        std::pow(a0 / fConfig->nuclearRadius, 2.) *
        (pqSurf[0] * pqSurf[0] + pqSurf[1] * pqSurf[1]);

      surfaceData.push_back(e * e0);
      surfaceData.push_back(probabilitySurface);

      if (!fConfig->writeWF)
        continue;
      std::string boundWFFileName = CreateFileName("radial", 0, iN, iK);
      // if (fConfig->outputDirectory == ""){
      //    boundWFFileName =
      //     "../Nuclei/" + fConfig->nucleusName + "/" + fConfig->nucleusName + "_" +
      //     fConfig->wfType + "_" + nameSeed + "_screening" +
      //     std::to_string(fConfig->applyScreening) + "_" + fConfig->processName +
      //     "_n" + std::to_string(iN) + "_k" + std::to_string(iK) + ".dat";
      // }
      // else{
      //   boundWFFileName = fConfig->outputDirectory + "/" + fConfig->nucleusName + "_" +
      //     fConfig->wfType + "_" + nameSeed + "_screening" +
      //     std::to_string(fConfig->applyScreening) + "_" + fConfig->processName +
      //     "_n" + std::to_string(iN) + "_k" + std::to_string(iK) + ".dat";
      // }
      
      int iMaxR = -1;
      for (int i = radwf.NGP - 1; i >= 0; i--) {
        if (abs(radwf.P[i]) > 1.E-35) {
          iMaxR = i;
          break;
        }
      }

      WriteRadialWFFile(iMaxR, e*e0, boundWFFileName.data());
    }
  }

  WriteSurfWFLine(0., surfaceData, surfaceWFFile);
}

void radial_interface::WriteRadialWFFile(int irMax, double energy, const std::string &boundWFFileName, double phase) {
  std::ofstream file(boundWFFileName.data(), std::ios::out | std::ios::binary);
  

  BoundRadialWFFileHeader hBound;
  ScatteringRadialWFFileHeader hScatter;
  if (fConfig->wfType.find("scattering") != std::string::npos){
    hScatter.energy = energy;
    hScatter.phase = phase;
    file.write((char*) &hScatter, sizeof(ScatteringRadialWFFileHeader));
  }
  else{
    hBound.energy = energy;
    file.write((char*) &hBound, sizeof(BoundRadialWFFileHeader));
  }

  RadialWFFileLine line;
  for (int iR = 0; iR < irMax; iR++){
    line.r = radwf.RAD[iR];
    line.p = radwf.P[iR];
    line.q = radwf.Q[iR];
    file.write((char*) &line, sizeof(RadialWFFileLine));
  }
  file.close();
}
