#include <cmath>

#include "functions.h"
#include "constants.h"

double diff_xsec(double ene_dm, double ene_r, double m_dm, double mA)
{
  const double overallFactor = 4.0 * M_PI * epsilon * epsilon * alpha_D * alpha / ( ene_dm * ene_dm - m_dm * m_dm );
  const double denominator = ( 2.0 * m_e * ene_r + mA * mA );
  return overallFactor * ( 2.0 * ene_dm * m_e * (ene_dm - ene_r) - m_dm * m_dm * ene_r ) / (denominator * denominator);
}

double myroot_diff_xsec(double* x, double* p)
{
  double mA = p[2];
  const double overallFactor = 4.0 * M_PI * epsilon * epsilon / ( p[0] * p[0] - p[1] * p[1] ) * alpha_D * alpha * (-1) * ( 2*m_e*p[0] + p[1]*p[1] ) / ( 4.*m_e*m_e );
  double diff = overallFactor * ( x[0] - ( 2.0*m_e*p[0]*p[0] / (2.0*m_e*p[0]+p[1]*p[1]) ) ) / ( ( x[0] + 0.5*mA * mA / m_e ) * ( x[0] + 0.5*mA * mA / m_e ) );
  return diff;
}

double xsec_analytic_integration(double ene_dm, double m_dm, double mA, double ene_i, double ene_f)
{
  double a = 2.0 * ene_dm * ene_dm * m_e / ( 2.0 * ene_dm * m_e + m_dm * m_dm);     // kinematic forbidden energy limit
  double b = 0.5 * ( mA * mA ) / m_e;
  double e_thr = ene_i;
  const double overallFactor = 4.0 * M_PI * epsilon * epsilon / ( ene_dm * ene_dm - m_dm * m_dm ) * alpha_D * alpha * (-1 * a);
  double integral = log( (a + b)/( e_thr + b) ) + (a+b)*( 1./(a + b) - 1./(e_thr + b ) );

  return overallFactor*integral;
}

double xsec_numerical_integration(double ene_dm, double m_dm, double mA, double ene_i, double ene_f)
{
  double e_0 = 2.0 * ene_dm * ene_dm * m_e / ( 2.0 * ene_dm * m_e + m_dm * m_dm );
  int n = 10000;
  double e_thr = ene_i;
  double h = (e_0 - e_thr)/n;
  double x;
  double param[3] = {ene_dm, m_dm, mA};
  double sum = 0.5 * ( myroot_diff_xsec(&e_thr, param) + myroot_diff_xsec(&e_0, param) );

  for(int i = 0; i < n; i++)
  {
    x = e_thr + i * h;

    sum += myroot_diff_xsec(&x, param);
  }
  return sum * h;
}
