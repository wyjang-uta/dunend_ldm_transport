#include <iostream>
#include <fstream>
#include "physconst.h"
#include "func.h"

/// input constants (data card)
/// physical variables
const double dm_ratio = 3.0;                            // m_A' = 3 * m_dm
/// scan range
const double mdm_minimum = 3.0;                         // mass scan window minimum in MeV
const double mdm_maximum = 300.0;                       // mass scan window maximum in MeV
const double e_thr = 30.0;                              // recoiled electron threshold energy in MeV

// temporary constants
const double energy_dm = 100.0;
const double energy_recoil = 50.0;
const double mass_dm = 1.0;
const double e_0 = 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm );

int main(int argc, char* argv[])
{
  double mA = dm_ratio * mass_dm;
  std::cout << "differential cross section with inputs: " << std::endl;
  std::cout << "Dark matter energy = " << energy_dm << " (MeV)" << std::endl;
  std::cout << "Dark matter mass = " << mass_dm << " (MeV)" << std::endl;
  std::cout << "Recoiled electron energy = " << energy_recoil << " (MeV)" << std::endl;
  std::cout << "Detector threshold = " << e_thr << " (MeV)" << std::endl;
  std::cout << "result (myfunc) = " << diff_xsec(energy_dm, energy_recoil, mass_dm, mA) << std::endl;
  //std::cout << "result (TF1) = " << fXsec->Eval(energy_recoil) << std::endl;

  double x = energy_recoil;
  double par[3] = {energy_dm, mass_dm, mA};
  std::cout << "result (myroot_diff_xsec) = " << myroot_diff_xsec(&x, par) << std::endl;
  double e_0 = 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm );
  std::cout << "The zero of the differential cross section is = " << e_0 << std::endl;
  std::cout << "Total cross section (analytic)= " << xsec_analytic_integration(energy_dm, mass_dm, mA, e_thr, e_0) << std::endl;
  std::cout << "Total cross section (numerical integration)= " << xsec_numerical_integration(energy_dm, mass_dm, mA, e_thr, e_0) << std::endl;

	return 0;
}

