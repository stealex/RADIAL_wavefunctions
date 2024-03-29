#ifndef MATHTOOLS_HH
#define MATHTOOLS_HH 1

#include <cmath>
#include <vector>

class math_tools {
public:
  static math_tools *GetInstance();
  void ComputeEnergyPoints(int nPoints, double minVal, double maxVal,
                           std::vector<double> &energyPoints);
  void GLNodesAndWeights(int nPoints);
  void Legangl(int n, double x, std::vector<double> &fx);

  double ChargedSpherePot(double z, double rSphere, double r);

  void ComputeScreeningFunction(int nPoints, double xMax, int nMaxCoeffs = 100,
                                bool dump = false);
  void ComputeExpansionCoefficients(std::vector<double> &coeffs);

  static std::vector<double> screenXvals;
  static std::vector<double> screenPhiVals;

  double InitialElectronDensity(double r, double b);
  double FermiProtonDensity(double r, double nuclearRadius, double diffusivityParam);

private:
  math_tools();
  virtual ~math_tools();
  static math_tools *fInstance;

  std::vector<double> glNodes;
  std::vector<double> glWeights;
};
#endif