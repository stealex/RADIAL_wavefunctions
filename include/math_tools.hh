#ifndef MATHTOOLS_HH
#define MATHTOOLS_HH 1

#include <vector>
class math_tools{
    public:
        static math_tools *GetInstance();
        void ComputeEnergyPoints(int nPoints, double maxVal, std::vector<double> &energyPoints);
        void GLNodesAndWeights(int nPoints);
        void Legangl(int n, double x, std::vector<double> &fx);

        double ChargedSpherePot(double z, double rSphere, double r);

        void ComputeScreeningFunction(int nPoints, double xMax, int nMaxCoeffs=100, bool dump=false);
        void ComputeExpansionCoefficients(std::vector<double> &coeffs);

        std::vector<double> screenXvals;
        std::vector<double> screenPhiVals;

    private:
        math_tools();
        virtual ~math_tools();
        static math_tools *fInstance;

        std::vector<double> glNodes;
        std::vector<double> glWeights;

};
#endif