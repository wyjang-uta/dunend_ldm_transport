#include <iostream>
#include <fstream>
#include <vector>

#include "Math/LorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"

#include "physconst.h"
#include "func.h"

/// input constants (data card)
/// physical variables
const double dm_ratio = 3.0;                            // m_A' = 3 * m_dm
/// scan range
const double mdm_minimum = 3.0;                         // mass scan window minimum in MeV
const double mdm_maximum = 300.0;                       // mass scan window maximum in MeV
const double e_thr = 30.0;                              // recoiled electron threshold energy in MeV
const char* pi0_infile = "dune_pi_zero_pot1e6_px_py_pz_E.txt";

// temporary constants
const double energy_dm = 100.0;                         // unit in MeV
const double energy_recoil = 50.0;                      // unit in MeV
const double mass_dm = 1.0;                             // unit in MeV
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

  // Loading pi0 file
  // - build 4-vector container

  std::vector<TLorentzVector> pi0s;
  TLorentzVector pi0;
  std::ifstream infile_pi0;
  const size_t bufferSize = 1024 * 1024; // size of the buffer (1 MB)
  infile_pi0.open(pi0_infile, std::ios::in | std::ios::binary);
  if( !infile_pi0.is_open() )
  {
    std::cerr << "Failed to open file: " << pi0_infile << '\n';
    return 1;
  }
  std::vector<char> buffer(bufferSize);
  infile_pi0.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
  double px, py, pz, E;
  long lineCount = 0;
  TH1D* hpizero = new TH1D("hpizero", "pi0 energy spectrum", 1200, 0, 120000);
  while( infile_pi0 >> px >> py >> pz >> E )
  {
    pi0.SetPxPyPzE(px, py, pz, E);
    pi0s.push_back(pi0);
    hpizero->Fill(pi0.E());
    if( lineCount % 100000 == 0 ) std::cout << '\r' << "Loading pi0 4-vector data from " << pi0_infile << ", " << lineCount << " particles have been loaded." << std::flush;
    lineCount++;
  }
  std::cout << '\r' << "Loading pi0 4-vector data from " << pi0_infile << ", " << lineCount << " particles have been successfully loaded." << std::endl;

  // Decaying pi0s' to obtain A'
  TGenPhaseSpace decay;
  TLorentzVector pi0_copy;
  TLorentzVector* darkPhoton;
  std::vector<TLorentzVector> darkPhotons;
  double masses[2] = {mA, 0}; // dark photon and massless SM photon
  long decayCount_pi0 = 0;
  TH1D* hAprime = new TH1D("hAprime", "dark photon energy spectrum", 1200, 0, 120000);
  // initialize random seed
  TRandom3 seed_pi0_decay(12345);
  for( auto& pi0_index : pi0s )
  {
    decay.SetDecay(pi0_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons.push_back(*darkPhoton);
    hAprime->Fill(darkPhoton->E());
    decayCount_pi0++;
  }
  std::cout << decayCount_pi0 << " neutral pions got decayed by TGenPhaseSpace code.\n";

  // Decaying A's to obtain phis
  masses[0] = mass_dm;
  masses[1] = mass_dm;
  long decayCount_phi = 0;
  // initialize random seed
  TRandom3 seed_darkPhoton_decay(1);
  std::vector<TLorentzVector> darkMatters;
  TLorentzVector* dm_particle1;
  TLorentzVector* dm_particle2;
  TH1D* hphi = new TH1D("hphi", "dark matter flux", 1200, 0, 120000);
  for( auto& darkPhoton_index : darkPhotons )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    //
    darkMatters.push_back(*dm_particle1);
    darkMatters.push_back(*dm_particle2);
    hphi->Fill(dm_particle1->E());
    hphi->Fill(dm_particle2->E());
    decayCount_phi++;
  }
  std::cout << decayCount_phi << " x 2 light dark mater particles created by TGenPhaseSpace code.\n";


  TFile fOutput("proto_output.root","RECREATE");
  hpizero->Write();
  hAprime->Write();
  hphi->Write();
  fOutput.Close();



  // Calculating expected number of signals
  //


  // Producing sensitivity curve
  //

  // Writing output
  //

	return 0;
}

