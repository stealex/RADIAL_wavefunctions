#ifndef PHYSICALCONSTANTS_HH
#define PHYSICALCONSTANTS_HH 1

namespace physical_constants {

double const pi = 3.1415926535897932;
double const a0 = 0.529177E5; // borh radius in fm
double const e0 = 27.2114E-6; // Hartree energy in MeV
double const electronMass = 0.510998950; // in MeV
double const hc = 197.32698950; // MeV*fm
double const fineStructure = 1. / 137.036; 
double const normalUnits = hc / (a0 * e0);
double const thomasFermiRadius = 0.88534; // prefactor

} // namespace physical_constants

#endif
