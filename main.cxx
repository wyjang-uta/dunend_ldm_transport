#include <iostream>
#include <fstream>
#include <vector>

#include "Math/LorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"

#include "constants.h"
#include "func.h"
#include "Vector3D.h"
#include "RayIntersect.h"

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
  //std::cout << "differential cross section with inputs: " << std::endl;
  //std::cout << "Dark matter energy = " << energy_dm << " (MeV)" << std::endl;
  //std::cout << "Dark matter mass = " << mass_dm << " (MeV)" << std::endl;
  //std::cout << "Recoiled electron energy = " << energy_recoil << " (MeV)" << std::endl;
  //std::cout << "Detector threshold = " << e_thr << " (MeV)" << std::endl;
  //std::cout << "result (myfunc) = " << diff_xsec(energy_dm, energy_recoil, mass_dm, mA) << std::endl;
  //std::cout << "result (TF1) = " << fXsec->Eval(energy_recoil) << std::endl;

  //double x = energy_recoil;
  //double par[3] = {energy_dm, mass_dm, mA};
  //std::cout << "result (myroot_diff_xsec) = " << myroot_diff_xsec(&x, par) << std::endl;
  //double e_0 = 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm );
  //std::cout << "The zero of the differential cross section is = " << e_0 << std::endl;
  //std::cout << "Total cross section (analytic)= " << xsec_analytic_integration(energy_dm, mass_dm, mA, e_thr, e_0) << std::endl;
  //std::cout << "Total cross section (numerical integration)= " << xsec_numerical_integration(energy_dm, mass_dm, mA, e_thr, e_0) << std::endl;

  // Load Geometry
  Vector3D v0(nd_location_vertex[0][0], nd_location_vertex[0][1], nd_location_vertex[0][2]);
  Vector3D v1(nd_location_vertex[1][0], nd_location_vertex[1][1], nd_location_vertex[1][2]);
  Vector3D v2(nd_location_vertex[2][0], nd_location_vertex[2][1], nd_location_vertex[2][2]);
  Vector3D v3(nd_location_vertex[3][0], nd_location_vertex[3][1], nd_location_vertex[3][2]);
  Vector3D v4(nd_location_vertex[4][0], nd_location_vertex[4][1], nd_location_vertex[4][2]);
  Vector3D v5(nd_location_vertex[5][0], nd_location_vertex[5][1], nd_location_vertex[5][2]);
  Vector3D v6(nd_location_vertex[6][0], nd_location_vertex[6][1], nd_location_vertex[6][2]);
  Vector3D v7(nd_location_vertex[7][0], nd_location_vertex[7][1], nd_location_vertex[7][2]);
  Box box(v0, v1, v2, v3, v4, v5, v6, v7);

  // Print Geometry
  std::cout << "nd_vertex 1: (" << nd_location_vertex[0][0] << ", " << nd_location_vertex[0][1] << ", " << nd_location_vertex[0][2] << ").\n";
  std::cout << "nd_vertex 2: (" << nd_location_vertex[1][0] << ", " << nd_location_vertex[1][1] << ", " << nd_location_vertex[1][2] << ").\n";
  std::cout << "nd_vertex 3: (" << nd_location_vertex[2][0] << ", " << nd_location_vertex[2][1] << ", " << nd_location_vertex[2][2] << ").\n";
  std::cout << "nd_vertex 4: (" << nd_location_vertex[3][0] << ", " << nd_location_vertex[3][1] << ", " << nd_location_vertex[3][2] << ").\n";
  std::cout << "nd_vertex 5: (" << nd_location_vertex[4][0] << ", " << nd_location_vertex[4][1] << ", " << nd_location_vertex[4][2] << ").\n";
  std::cout << "nd_vertex 6: (" << nd_location_vertex[5][0] << ", " << nd_location_vertex[5][1] << ", " << nd_location_vertex[5][2] << ").\n";
  std::cout << "nd_vertex 7: (" << nd_location_vertex[6][0] << ", " << nd_location_vertex[6][1] << ", " << nd_location_vertex[6][2] << ").\n";
  std::cout << "nd_vertex 8: (" << nd_location_vertex[7][0] << ", " << nd_location_vertex[7][1] << ", " << nd_location_vertex[7][2] << ").\n";

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
  TH1D* hpizero = new TH1D("hpizero", "#pi_{0} energy spectrum;#pi_{0} Energy (MeV);Entries", nbins_E, fEx_i, fEx_f);
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
  double mA = dm_ratio * mass_dm;
  double masses[2] = {mA, 0}; // dark photon and massless SM photon
  long decayCount_pi0 = 0;
  TH1D* hAprime = new TH1D("hAprime", "Dark photon energy spectrum;A' Energy (MeV);Entries", nbins_E, fEx_i, fEx_f);
  // initialize random seed
  TRandom3 seed_pi0_decay(12345);
  for( auto& pi0_index : pi0s )
  {
    decay.SetDecay(pi0_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons.push_back(*darkPhoton);
    hAprime->Fill(darkPhoton->E());
    if( decayCount_pi0 % 100000 == 0 ) std::cout << '\r' << "Decaying pi0s' into photon + dark photon pair ... " << decayCount_pi0 << " particles have been loaded. (" << 100*(double)decayCount_pi0/lineCount << ")" << std::flush;
    decayCount_pi0++;
  }
  std::cout << '\r' << "Decaying pi0s' into photon + dark photon pair ... DONE. " << decayCount_pi0 << " pions decayed." << std::endl;

  // Decaying A's to obtain phis
  masses[0] = mass_dm;
  masses[1] = mass_dm;
  long decayCount_phi = 0;
  long decayCount_phi_accepted = 0;
  // initialize random seed
  TRandom3 seed_darkPhoton_decay(1);
  std::vector<TLorentzVector> darkMatters;
  TLorentzVector* dm_particle1;
  TLorentzVector* dm_particle2;
  TH1D* hphi = new TH1D("hphi", "Dark matter energy spectrum (total);#phi Energy (MeV);Entries", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_acc = new TH1D("hphi_acc", "Dark matter energy spectrum (accepted);#phi Energy (MeV);Entries", nbins_E, fEx_i, fEx_f);
  Vector3D rayO(0.0, 0.0, 0.0);
  Vector3D rayV(0.0, 0.0, 0.0);
  bool interCube;
  for( auto& darkPhoton_index : darkPhotons )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    rayV.x = dm_particle1->Px();
    rayV.y = dm_particle1->Py();
    rayV.z = dm_particle1->Pz();
    hphi->Fill(dm_particle1->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle1);
      hphi_acc->Fill(dm_particle1->E());
      decayCount_phi_accepted++;
    }
    rayV.x = dm_particle2->Px();
    rayV.y = dm_particle2->Py();
    rayV.z = dm_particle2->Pz();
    hphi->Fill(dm_particle2->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle2);
      hphi_acc->Fill(dm_particle2->E());
      decayCount_phi_accepted++;
    }
    if( decayCount_phi % 100000 == 0 ) std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been loaded. (" << 100*(double)decayCount_phi/lineCount << ")" << std::flush;
    //
    decayCount_phi++;
  }
  std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been decayed. (" << 100*(double)decayCount_phi/lineCount << ")" << std::endl;
  std::cout << decayCount_phi_accepted << " light dark mater particles accepted within DUNE ND fiducial volume.\n";
  std::cout << "Geometrical acceptance : " << (double)decayCount_phi_accepted / decayCount_phi * 4.0 * M_PI * dist_mc0_nd * dist_mc0_nd << "m^2" << '\n';

  // Calculate expected number of signals
  //


  // Producing sensitivity curve
  //

  // Writing output
  //
  TFile fOutput("proto_output.root","RECREATE");
  hpizero->Write();
  hAprime->Write();
  hAprime->Write();
  hphi->Write();
  hphi_acc->Write();
  TH1D* hphi_eff = (TH1D*)hphi_acc->Clone("hphi_eff");
  hphi_eff->Divide(hphi);
  hphi_eff->SetTitle("Efficiency");
  hphi_eff->GetYaxis()->SetTitle("Efficiency");
  hphi_eff->Write();

  fOutput.Close();

	return 0;
}

