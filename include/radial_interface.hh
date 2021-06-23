#ifndef RADIALINTERFACE_HH
#define RADIALINTERFACE_HH 1

#include <string>
#include <vector>
#include <fstream>
#include "config_utility.hh"

class radial_interface{
    public:
        static radial_interface *GetInstance();

        void Initialize();
        void ObtainPotential(double zPotential);
        void ApplyScreening(int nPoints = 1000);

        void SolveDirac();
        void SolveScatteringStates();
        void SolveBoundStates();

        void WriteSurfWFLine(double e, std::vector<double> vals, std::ofstream &file);
        void WriteRadialWFFile(int irMax,std::ofstream &file);
        void FindPQSurface(std::vector<double> &pq);

    private:
        radial_interface();
        virtual ~radial_interface();
        static radial_interface *fInstance;
        config_utility *fConfig;

        std::vector<double> energyPoints;
        std::vector<double> rValues;
        std::vector<double> potValues;
        std::vector<double> rvValues;
};
#endif