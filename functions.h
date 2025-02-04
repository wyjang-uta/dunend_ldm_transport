#ifndef FUNC_H
#define FUNC_H

double diff_xsec(double ene_dm, double ene_r, double m_dm, double mA);
double myroot_diff_xsec(double* x, double* param);
double xsec_analytic_integration(double ene_dm, double m_dm, double mA, double ene_i, double ene_f);
double xsec_numerical_integration(double ene_dm, double m_dm, double mA, double ene_i, double ene_f);

#endif

