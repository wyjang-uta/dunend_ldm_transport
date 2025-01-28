#include <iostream>
#include "physconst.h"

/// input constants
/// physical variables
const double dm_ratio = 3.0;                            // m_A' = 3 * m_dm
const double energy_dm = 100.0;
const double energy_recoil = 50.0;
const double mass_dm = 1.0;
const double e_0 = 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm );
/// scan range
const double mdm_minimum = 3.0;                         // mass scan window minimum in MeV
const double mdm_maximum = 300.0;                       // mass scan window maximum in MeV
const double e_thr = 30.0;                              // recoiled electron threshold energy in MeV

double myroot_diff_xsec(double* x, double* param);
double diff_xsec(double ene_dm, double ene_r, double m_dm);
double xsec_analytic_integration(double ene_dm, double m_dm);
double xsec_numerical_integration(double ene_dm, double m_dm);

int main()
{
  TF1* fXsec = new TF1("fXsec", myroot_diff_xsec, e_thr, e_0, 2);
  fXsec->GetXaxis()->SetTitle("Recoiled Electron Energy (MeV)");
  fXsec->GetYaxis()->SetTitle("Differential Cross-Section (MeV^{-2})");
  fXsec->GetYaxis()->SetLabelFont(22);
  fXsec->GetYaxis()->SetLabelSize(0.05);
  fXsec->GetYaxis()->SetTitleFont(22);
  fXsec->GetYaxis()->SetTitleSize(0.05);
  fXsec->SetParameters(energy_dm, mass_dm);
  fXsec->SetNpx(10000);


  std::cout << "differential cross section with inputs: " << std::endl;
  std::cout << "Dark matter energy = " << energy_dm << " (MeV)" << std::endl;
  std::cout << "Dark matter mass = " << mass_dm << " (MeV)" << std::endl;
  std::cout << "Recoiled electron energy = " << energy_recoil << " (MeV)" << std::endl;
  std::cout << "Detector threshold = " << e_thr << " (MeV)" << std::endl;
  std::cout << "result (myfunc) = " << diff_xsec(energy_dm, energy_recoil, mass_dm) << std::endl;
  std::cout << "result (TF1) = " << fXsec->Eval(energy_recoil) << std::endl;

  double x = energy_recoil;
  double par[2] = {energy_dm, mass_dm};
  std::cout << "result (myroot_diff_xsec) = " << myroot_diff_xsec(&x, par) << std::endl;
  std::cout << "The zero of the differential cross section is = " << 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm ) << std::endl;
  std::cout << "Total cross section (analytic)= " << xsec_analytic_integration(energy_dm, mass_dm) << std::endl;
  std::cout << "Total cross section (numerical integration)= " << xsec_numerical_integration(energy_dm, mass_dm) << std::endl;
  std::cout  << "Total cross section (TF1) = " << fXsec->Integral(e_thr, e_0) << std::endl;

	return 0;
}

double diff_xsec(double ene_dm, double ene_r, double m_dm)
{
  const double overallFactor = four_pi_epsilon_squared * alpha_D * alpha / ( ene_dm * ene_dm - m_dm * m_dm );
  const double denominator = ( 2.0 * m_e * ene_r + dm_ratio * dm_ratio * m_dm * m_dm );
  return overallFactor * ( 2.0 * ene_dm * m_e * (ene_dm - ene_r) - m_dm * m_dm * ene_r ) / (denominator * denominator);
}

double xsec_analytic_integration(double ene_dm, double m_dm)
{
  double mA = dm_ratio * m_dm;
  double a = e_0;
  double b = 0.5 * ( mA * mA ) / m_e;
  const double overallFactor = 4.0 * M_PI * epsilon * epsilon / ( ene_dm * ene_dm - m_dm * m_dm ) * alpha_D * alpha * (-1 * a);
  double integral = log( (e_0 + b)/( e_thr + b) ) + (a+b)*( 1./(e_0 + b) - 1./(e_thr + b ) );

  return overallFactor*integral;
}

double myroot_diff_xsec(double* x, double* p)
{
  double mA = dm_ratio * p[1];
  const double overallFactor = 4.0 * M_PI * epsilon * epsilon / ( p[0] * p[0] - p[1] * p[1] ) * alpha_D * alpha * (-1) * (2*m_e*p[0] + p[1]*p[1]) / (4.*m_e*m_e);
  double diff = overallFactor * ( x[0] - (2.0*m_e*p[0]*p[0] / (2.0*m_e*p[0]+p[1]*p[1])) )/(( x[0] + 0.5*mA * mA / m_e)*( x[0] + 0.5*mA * mA / m_e));
  return diff;
}

double xsec_numerical_integration(double ene_dm, double m_dm)
{
  double sum = 0.0;
  double e_0 = 2.0 * ene_dm * ene_dm * m_e / ( 2.0 * ene_dm * m_e + m_dm * m_dm );
  int n = 10000;
  double h = (e_0 - e_thr)/n;
  double x;
  double param[2] = {ene_dm, m_dm};

  for(int i = 0; i < n; i++)
  {
    x = e_thr + i * h;

    sum += myroot_diff_xsec(&x, param) * h;
  }
  return sum;
}
