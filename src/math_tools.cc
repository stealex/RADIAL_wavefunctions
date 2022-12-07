#include "math_tools.hh"
#include "physical_constants.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

using namespace physical_constants;

math_tools *math_tools::fInstance = 0;
std::vector<double> math_tools::screenXvals = {};
std::vector<double> math_tools::screenPhiVals = {};

math_tools *math_tools::GetInstance() {
  if (!fInstance)
    fInstance = new math_tools();
  return fInstance;
}

math_tools::math_tools() : glNodes(0), glWeights(0) {}

math_tools::~math_tools() {
  glNodes.clear();
  glWeights.clear();
}

void math_tools::ComputeEnergyPoints(int nPoints, double minValue,
                                     double maxValue,
                                     std::vector<double> &energyPoints) {
  if (energyPoints.size() != 0)
    energyPoints.clear();
  energyPoints.resize(nPoints);

  double eMaxLog = std::log10(maxValue);
  double eMinLog = std::log10(minValue);
  double de = (eMaxLog - eMinLog) / (nPoints - 1.);

  for (int i = 0; i < nPoints; i++) {
    double elog = eMinLog + i * de;
    energyPoints[i] = std::pow(10., elog);
  }
}

void math_tools::GLNodesAndWeights(int nPoints) {
  std::vector<double> fx(nPoints + 1);
  std::vector<double> f1x(nPoints);
  double fpx;

  glNodes.clear();
  glWeights.clear();

  glNodes.resize(nPoints);
  glWeights.resize(nPoints);
  for (int i = 0; i < nPoints; i++) {
    glNodes[i] = std::cos(pi * (i + 1.0 - 0.25) / (nPoints + 0.5));
    for (int j = 0; j < 20; j++) {
      Legangl(nPoints + 1, glNodes[i], fx);
      Legangl(nPoints, glNodes[i], f1x);
      fpx = nPoints * (glNodes[i] * fx[nPoints] - f1x[nPoints - 1]) /
            (std::pow(glNodes[i], 2) - 1.);
      glNodes[i] = glNodes[i] - fx[nPoints] / fpx;
      if (abs(fx[nPoints]) < 1E-12)
        break;
    }
    glWeights[i] = 2. / ((1. - glNodes[i] * glNodes[i]) * fpx * fpx);
  }
}

void math_tools::Legangl(int n, double x, std::vector<double> &fx) {
  fx[0] = 1.;
  fx[1] = x;
  for (int i = 1; i < n; i++) {
    fx[i + 1] = ((i * 2. + 1.) * x * fx[i] - i * fx[i - 1]) / (i + 1.);
  }
}

double math_tools::ChargedSpherePot(double z, double rSphere, double r) {
  if (r <= rSphere) {
    return z * fineStructure * hc / (2. * rSphere) *
           (3. - std::pow(r / rSphere, 2.));
  } else {
    return z * fineStructure * hc / r;
  }

  return 0.;
}

void math_tools::ComputeScreeningFunction(int nPoints, double xMax,
                                          int nMaxCoeff, bool dump) {
  std::vector<double> expansCoeffs(nMaxCoeff);
  ComputeExpansionCoefficients(expansCoeffs);

  std::vector<double> xGL(31);
  std::vector<double> wGL(31);
  GLNodesAndWeights(31);

  for (int i = 0; i < 31; i++) {
    xGL[i] = glNodes[31 - i - 1];
    wGL[i] = glWeights[31 - i - 1];
  }

  std::ofstream screeningFile;
  if (dump)
    screeningFile.open("screening_test.dat", std::ofstream::trunc);

  bool isLast = 0;
  for (int i = 0; i < nPoints; i++) {
    double tMax = i / (nPoints - 1.);

    double sum = 0.;
    double tVals[31];
    double uVals[31];
    double wVals[31];
    for (int j = 0; j < 31; j++) {
      tVals[j] = 0.5 * tMax * xGL[j] + 0.5 * tMax;

      uVals[j] = 0.0;
      for (int m = 0; m < nMaxCoeff; m++) {
        uVals[j] += expansCoeffs[m] * std::pow(1. - tVals[j], m);
      }

      wVals[j] =
          -6. * uVals[j] * tVals[j] / (1. - tVals[j] * tVals[j] * uVals[j]);
      sum += 0.5 * tMax * wGL[j] * wVals[j];
    }

    double x = std::pow(144., 1. / 3.) * tMax * tMax * std::exp(-sum / 3.);
    if (x > xMax)
      isLast = true;

    double phi = std::exp(sum);
    if (phi < 1E-15)
      phi = 0.;

    screenXvals.push_back(x);
    screenPhiVals.push_back(phi);

    if (!dump)
      continue;
    screeningFile.precision(6);
    screeningFile.width(14);
    screeningFile << std::scientific << x;
    screeningFile.precision(6);
    screeningFile.width(14);
    screeningFile << std::scientific << phi;
    screeningFile << std::endl;

    if (isLast)
      break;
  }
  if (dump)
    screeningFile.close();
}

void math_tools::ComputeExpansionCoefficients(std::vector<double> &coeffs) {
  int nCoeffs = (int)coeffs.size();

  coeffs[0] = 1.0;
  coeffs[1] = 9.0 - std::sqrt(73.0);

  for (int m = 2; m < nCoeffs; m++) {
    coeffs[m] = 1.0 / (2.0 * (m + 8.0) - (m + 1.0) * coeffs[1]);

    double sum = 0.0;
    for (int n = 1; n < m - 2; n++) {
      sum += coeffs[m - n] *
             ((n + 1.) * coeffs[n + 1] - 2. * (n + 4.) * coeffs[n] +
              (n + 7.) * coeffs[n - 1]);
    }
    sum += coeffs[m - 1] * ((m + 7.) - 2 * (m + 3.) * coeffs[1]) +
           coeffs[m - 2] * ((m + 6.) * coeffs[1]);
    coeffs[m] = coeffs[m] * sum;
  }
}

double math_tools::InitialElectronDensity(double r, double b) {
  return (1. / 0.07957747154594769) / (4 * pi * r * b * b) *
         (3.60 * std::exp(-6. * r / b) + 0.792 * std::exp(-1.2 * r / b) +
          0.0315 * std::exp(-0.3 * r / b));
}

double math_tools::FermiProtonDensity(double r, double nuclearRadius, double diffusivityParam = 0.546){
  return 1./(std::exp((r-nuclearRadius)/diffusivityParam) + 1.);
}