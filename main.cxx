#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <Math/LorentzVector.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>

#include <G4GDMLParser.hh>

#include "constants.h"
#include "functions.h"
#include "vector3D.h"
#include "rayIntersect.h"
#include "twoBodyDecayCalculator.h"

/// input constants (data card)
/// physical variables
const double dm_ratio = 3.0;                            // m_A' = 3 * m_dm
/// scan range
const double mdm_minimum = 3.0;                         // mass scan window minimum in MeV
const double mdm_maximum = 300.0;                       // mass scan window maximum in MeV
const double e_thr = 30.0;                              // recoiled electron threshold energy in MeV
const char* pi0_target_infile = "dune_pi_zero_target_pot1e6_px_py_pz_E.txt";
const char* eta_target_infile = "dune_eta_target_pot1e6_px_py_pz_E.txt";
const char* pi0_dump_infile = "dune_pi_zero_dump_pot1e4_px_py_pz_E.txt";
const char* eta_dump_infile = "dune_eta_dump_pot1e4_px_py_pz_E.txt";

// temporary constants
const double energy_dm = 100.0;                         // unit in MeV
const double energy_recoil = 50.0;                      // unit in MeV
const double mass_dm = 1.0;                             // unit in MeV
const double e_0 = 2.0 * energy_dm * energy_dm * m_e / ( 2.0 * energy_dm * m_e + mass_dm * mass_dm );

int main(int argc, char* argv[])
{
  double px, py, pz, E;
  double mA = dm_ratio * mass_dm;

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
  vector3D v0(nd_location_vertex[0][0], nd_location_vertex[0][1], nd_location_vertex[0][2]);
  vector3D v1(nd_location_vertex[1][0], nd_location_vertex[1][1], nd_location_vertex[1][2]);
  vector3D v2(nd_location_vertex[2][0], nd_location_vertex[2][1], nd_location_vertex[2][2]);
  vector3D v3(nd_location_vertex[3][0], nd_location_vertex[3][1], nd_location_vertex[3][2]);
  vector3D v4(nd_location_vertex[4][0], nd_location_vertex[4][1], nd_location_vertex[4][2]);
  vector3D v5(nd_location_vertex[5][0], nd_location_vertex[5][1], nd_location_vertex[5][2]);
  vector3D v6(nd_location_vertex[6][0], nd_location_vertex[6][1], nd_location_vertex[6][2]);
  vector3D v7(nd_location_vertex[7][0], nd_location_vertex[7][1], nd_location_vertex[7][2]);
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

  // Loading meson files
  // - build 4-vector container

  /*
  std::vector<TLorentzVector> pi0s_target;
  std::vector<TLorentzVector> pi0s_dump;
  std::vector<TLorentzVector> etas_target;
  std::vector<TLorentzVector> etas_dump;
  TLorentzVector pi0;
  TLorentzVector eta;
  std::ifstream infile_pi0_target;
  std::ifstream infile_pi0_dump;
  std::ifstream infile_eta_target;
  std::ifstream infile_eta_dump;
  const size_t bufferSize = 1024 * 1024; // size of the buffer (1 MB)
  infile_pi0_target.open(pi0_target_infile, std::ios::in | std::ios::binary);
  if( !infile_pi0_target.is_open() )
  {
    std::cerr << "Failed to open file: " << pi0_target_infile << '\n';
    return 1;
  }
  infile_pi0_dump.open(pi0_dump_infile, std::ios::in | std::ios::binary);
  if( !infile_pi0_dump.is_open() )
  {
    std::cerr << "Failed to open file: " << pi0_dump_infile << '\n';
    return 1;
  }
  infile_eta_target.open(eta_target_infile, std::ios::in | std::ios::binary);
  if( !infile_eta_target.is_open() )
  {
    std::cerr << "Failed to open file: " << eta_target_infile << '\n';
    return 1;
  }
  infile_eta_dump.open(eta_dump_infile, std::ios::in | std::ios::binary);
  if( !infile_eta_dump.is_open() )
  {
    std::cerr << "Failed to open file: " << eta_dump_infile << '\n';
    return 1;
  }
  std::vector<char> buffer(bufferSize);
  infile_pi0_target.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
  infile_pi0_dump.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
  infile_eta_target.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
  infile_eta_dump.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
  long lineCount = 0;
  TH1D* hpizero_target = new TH1D("hpizero_target", ";#pi_{0} Energy (MeV);#pi_{0} / year /#Delta E", nbins_E, fEx_i, fEx_f);
  TH1D* hpizero_dump = new TH1D("hpizero_dump", ";#pi_{0} Energy (MeV);#pi_{0} / year /#Delta E", nbins_E, fEx_i, fEx_f);
  TH1D* heta_target = new TH1D("heta_target", ";#eta Energy (MeV);#eta / year /#Delta E", nbins_E, fEx_i, fEx_f);
  TH1D* heta_dump = new TH1D("heta_dump", ";#eta Energy (MeV);#eta / year /#Delta E", nbins_E, fEx_i, fEx_f);
  */

  twoBodyDecayCalculator dcal_pi0_target(std::string(pi0_target_infile), m_pi0, mA, 0.0);
  dcal_pi0_target.load();
  dcal_pi0_target.decay();

  twoBodyDecayCalculator dcal_pi0_dump(std::string(pi0_dump_infile), m_pi0, mA, 0.0);
  dcal_pi0_dump.load();
  dcal_pi0_dump.decay();

  twoBodyDecayCalculator dcal_eta_target(std::string(eta_target_infile), m_eta, mA, 0.0);
  dcal_eta_target.load();
  dcal_eta_target.decay();
  
  twoBodyDecayCalculator dcal_eta_dump(std::string(eta_dump_infile), m_eta, mA, 0.0);
  dcal_eta_dump.load();
  dcal_eta_dump.decay();

  const std::vector<TLorentzVector>& dark_photon_pi0_target = dcal_pi0_target.getDaughter1();
  const std::vector<TLorentzVector>& dark_photon_pi0_dump   = dcal_pi0_dump.getDaughter1();
  const std::vector<TLorentzVector>& dark_photon_eta_target = dcal_eta_target.getDaughter1();
  const std::vector<TLorentzVector>& dark_photon_eta_dump   = dcal_eta_dump.getDaughter1();

  twoBodyDecayCalculator dcal_pi0_dphoton_target(dark_photon_pi0_target, mA, mass_dm, mass_dm);
  twoBodyDecayCalculator dcal_pi0_dphoton_dump(dark_photon_pi0_dump, mA, mass_dm, mass_dm);
  twoBodyDecayCalculator dcal_eta_dphoton_target(dark_photon_eta_target, mA, mass_dm, mass_dm);
  twoBodyDecayCalculator dcal_eta_dphoton_dump(dark_photon_eta_dump, mA, mass_dm, mass_dm);

  dcal_pi0_dphoton_target.decay();
  dcal_pi0_dphoton_dump.decay();
  dcal_eta_dphoton_target.decay();
  dcal_eta_dphoton_dump.decay();

  // now export the dark matter spectra

  /*
  while( infile_pi0_target >> px >> py >> pz >> E )
  {
    pi0.SetPxPyPzE(px, py, pz, E);
    pi0s_target.push_back(pi0);
    hpizero_target->Fill(pi0.E());
    if( lineCount % 100000 == 0 ) std::cout << '\r' << "Loading pi0 4-vector data from " << pi0_target_infile << ", " << lineCount << " particles have been loaded." << std::flush;
    lineCount++;
  }
  std::cout << '\r' << "Loading pi0 from target 4-vector data from " << pi0_target_infile << ", " << lineCount << " particles have been successfully loaded." << std::endl;

  lineCount = 0;
  while( infile_pi0_dump >> px >> py >> pz >> E )
  {
    pi0.SetPxPyPzE(px, py, pz, E);
    pi0s_dump.push_back(pi0);
    hpizero_dump->Fill(pi0.E());
    if( lineCount % 100000 == 0 ) std::cout << '\r' << "Loading pi0 4-vector data from " << pi0_dump_infile << ", " << lineCount << " particles have been loaded." << std::flush;
    lineCount++;
  }
  std::cout << '\r' << "Loading pi0 from dump proton 4-vector data from " << pi0_dump_infile << ", " << lineCount << " particles have been successfully loaded." << std::endl;

  lineCount = 0;
  while( infile_eta_target >> px >> py >> pz >> E )
  {
    eta.SetPxPyPzE(px, py, pz, E);
    etas_target.push_back(pi0);
    heta_target->Fill(eta.E());
    if( lineCount % 100000 == 0 ) std::cout << '\r' << "Loading eta 4-vector data from " << eta_target_infile << ", " << lineCount << " particles have been loaded." << std::flush;
    lineCount++;
  }
  std::cout << '\r' << "Loading eta from target 4-vector data from " << pi0_target_infile << ", " << lineCount << " particles have been successfully loaded." << std::endl;

  lineCount = 0;
  while( infile_eta_dump >> px >> py >> pz >> e )
  {
    eta.setpxpypze(px, py, pz, e);
    etas_dump.push_back(pi0);
    heta_dump->fill(eta.e());
    if( linecount % 100000 == 0 ) std::cout << '\r' << "loading eta 4-vector data from " << eta_dump_infile << ", " << linecount << " particles have been loaded." << std::flush;
    linecount++;
  }
  std::cout << '\r' << "loading eta from dump proton 4-vector data from " << pi0_dump_infile << ", " << linecount << " particles have been successfully loaded." << std::endl;

  // Decaying pi0s_target' to obtain A'
  TGenPhaseSpace decay;
  TLorentzVector* darkPhoton;
  std::vector<TLorentzVector> darkPhotons_pi0_target;
  std::vector<TLorentzVector> darkPhotons_pi0_dump;
  std::vector<TLorentzVector> darkPhotons_eta_target;
  std::vector<TLorentzVector> darkPhotons_eta_dump;
  double masses[2] = {mA, 0}; // dark photon and massless SM photon
  long decayCount_pi0_target = 0;
  long decayCount_pi0_dump = 0;
  long decayCount_eta_target = 0;
  long decayCount_eta_dump = 0;
  TH1D* hAprime_pi0_target  = new TH1D("hAprime_pi0_target", "Dark photon energy spectrum;A' Energy (MeV);A' production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hAprime_pi0_dump = new TH1D("hAprime_pi0_dump", "Dark photon energy spectrum;A' Energy (MeV);A' production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hAprime_eta_target  = new TH1D("hAprime_eta_target", "Dark photon energy spectrum;A' Energy (MeV);A' production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hAprime_eta_dump = new TH1D("hAprime_eta_dump", "Dark photon energy spectrum;A' Energy (MeV);A' production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);

  // initialize random seed
  TRandom3 seed_pi0_decay(12345);
  for( auto& pi0_index : pi0s_target )
  {
    decay.SetDecay(pi0_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons_pi0_target.push_back(*darkPhoton);
    hAprime_pi0_target->Fill(darkPhoton->E());
    if( decayCount_pi0_target % 100000 == 0 ) std::cout << '\r' << "Decaying pi0s_target' into photon + dark photon pair ... " << decayCount_pi0_target << " particles have been loaded. (" << 100*(double)decayCount_pi0_target/lineCount << ")" << std::flush;
    decayCount_pi0_target++;
  }
  std::cout << '\r' << "Decaying pi0s_target' into photon + dark photon pair ... DONE. " << decayCount_pi0_target << " pions decayed." << std::endl;

  for( auto& pi0_index : pi0s_dump )
  {
    decay.SetDecay(pi0_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons_pi0_dump.push_back(*darkPhoton);
    hAprime_pi0_dump->Fill(darkPhoton->E());
    if( decayCount_pi0_dump % 100000 == 0 ) std::cout << '\r' << "Decaying pi0s_dump' into photon + dark photon pair ... " << decayCount_pi0_dump << " particles have been loaded. (" << 100*(double)decayCount_pi0_dump/lineCount << ")" << std::flush;
    decayCount_pi0_dump++;
  }
  std::cout << '\r' << "Decaying pi0s_dump' into photon + dark photon pair ... DONE. " << decayCount_pi0_dump << " pions decayed." << std::flush << std::endl;

  for( auto& eta_index : etas_target )
  {
    decay.SetDecay(eta_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons_eta_target.push_back(*darkPhoton);
    hAprime_eta_target->Fill(darkPhoton->E());
    if( decayCount_eta_target % 100000 == 0 ) std::cout << '\r' << "Decaying etas_target' into photon + dark photon pair ... " << decayCount_eta_target << " particles have been loaded. (" << 100*(double)decayCount_eta_target/lineCount << ")" << std::flush;
    decayCount_eta_target++;
  }
  std::cout << '\r' << "Decaying etas_target' into photon + dark photon pair ... DONE. " << decayCount_eta_target << " pions decayed." << std::endl;


  for( auto& eta_index : etas_dump )
  {
    decay.SetDecay(eta_index, 2, masses, "");
    decay.Generate();
    darkPhoton = decay.GetDecay(0);
    darkPhotons_eta_dump.push_back(*darkPhoton);
    hAprime_eta_dump->Fill(darkPhoton->E());
    if( decayCount_eta_dump % 100000 == 0 ) std::cout << '\r' << "Decaying etas_dump' into photon + dark photon pair ... " << decayCount_eta_dump << " particles have been loaded. (" << 100*(double)decayCount_eta_dump/lineCount << ")" << std::flush;
    decayCount_eta_dump++;
  }
  std::cout << '\r' << "Decaying etas_dump' into photon + dark photon pair ... DONE. " << decayCount_eta_dump << " pions decayed." << std::endl;

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
  TH1D* hphi_pi0_target      = new TH1D("hphi_pi0_target", "Dark matter energy spectrum (total);#phi Energy (MeV);#phi production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_pi0_target_acc  = new TH1D("hphi_pi0_target_acc", "Dark matter energy spectrum (accepted);#phi Energy (MeV);#phi in ND / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_pi0_dump     = new TH1D("hphi_pi0_dump", "Dark matter energy spectrum (total);#phi Energy (MeV);#phi production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_pi0_dump_acc = new TH1D("hphi_pi0_dump_acc", "Dark matter energy spectrum (accepted);#phi Energy (MeV);#phi in ND / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_eta_target      = new TH1D("hphi_eta_target", "Dark matter energy spectrum (total);#phi Energy (MeV);#phi production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_eta_target_acc  = new TH1D("hphi_eta_target_acc", "Dark matter energy spectrum (accepted);#phi Energy (MeV);#phi in ND / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_eta_dump     = new TH1D("hphi_eta_dump", "Dark matter energy spectrum (total);#phi Energy (MeV);#phi production / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  TH1D* hphi_eta_dump_acc = new TH1D("hphi_eta_dump_acc", "Dark matter energy spectrum (accepted);#phi Energy (MeV);#phi in ND / year /#Delta E (#epsilon = 1)", nbins_E, fEx_i, fEx_f);
  vector3D rayO(0.0, 0.0, 0.0);
  vector3D rayV(0.0, 0.0, 0.0);
  bool interCube;
  for( auto& darkPhoton_index : darkPhotons_pi0_target )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    rayV.setX(dm_particle1->Px());
    rayV.setY(dm_particle1->Py());
    rayV.setZ(dm_particle1->Pz());
    hphi_pi0_target->Fill(dm_particle1->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle1);
      hphi_pi0_target_acc->Fill(dm_particle1->E());
      decayCount_phi_accepted++;
    }
    rayV.setX(dm_particle2->Px());
    rayV.setY(dm_particle2->Py());
    rayV.setZ(dm_particle2->Pz());
    hphi_pi0_target->Fill(dm_particle2->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle2);
      hphi_pi0_target_acc->Fill(dm_particle2->E());
      decayCount_phi_accepted++;
    }
    if( decayCount_phi % 100000 == 0 ) std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been loaded. (" << 100*(double)decayCount_phi/lineCount << ")" << std::flush;
    //
    decayCount_phi++;
  }
  std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been decayed. (" << 100.*(double)decayCount_phi/lineCount << ")" << std::endl;
  std::cout << decayCount_phi_accepted << " light dark mater particles accepted within DUNE ND fiducial volume.\n";
  std::cout << "Geometrical acceptance : " << (double)decayCount_phi_accepted / decayCount_phi * 4.0 * M_PI * dist_mc0_nd * dist_mc0_nd << "m^2" << '\n';

  for( auto& darkPhoton_index : darkPhotons_pi0_dump )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    rayV.setX(dm_particle1->Px());
    rayV.setY(dm_particle1->Py());
    rayV.setZ(dm_particle1->Pz());
    hphi_pi0_dump->Fill(dm_particle1->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle1);
      hphi_pi0_dump_acc->Fill(dm_particle1->E());
      decayCount_phi_accepted++;
    }
    rayV.setX(dm_particle2->Px());
    rayV.setY(dm_particle2->Py());
    rayV.setZ(dm_particle2->Pz());
    hphi_pi0_dump->Fill(dm_particle2->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle2);
      hphi_pi0_dump_acc->Fill(dm_particle2->E());
      decayCount_phi_accepted++;
    }
    if( decayCount_phi % 100000 == 0 ) std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been loaded. (" << 100*(double)decayCount_phi/lineCount << ")" << std::flush;
    //
    decayCount_phi++;
  }
  std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been decayed. (" << 100.*(double)decayCount_phi/lineCount << ")" << std::endl;
  std::cout << decayCount_phi_accepted << " light dark mater particles accepted within DUNE ND fiducial volume.\n";
  std::cout << "Geometrical acceptance : " << (double)decayCount_phi_accepted / decayCount_phi * 4.0 * M_PI * dist_mc0_nd * dist_mc0_nd << "m^2" << '\n';

  for( auto& darkPhoton_index : darkPhotons_eta_target )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    rayV.setX(dm_particle1->Px());
    rayV.setY(dm_particle1->Py());
    rayV.setZ(dm_particle1->Pz());
    hphi_eta_target->Fill(dm_particle1->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle1);
      hphi_eta_target_acc->Fill(dm_particle1->E());
      decayCount_phi_accepted++;
    }
    rayV.setX(dm_particle2->Px());
    rayV.setY(dm_particle2->Py());
    rayV.setZ(dm_particle2->Pz());
    hphi_eta_target->Fill(dm_particle2->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle2);
      hphi_eta_target_acc->Fill(dm_particle2->E());
      decayCount_phi_accepted++;
    }
    if( decayCount_phi % 100000 == 0 ) std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been loaded. (" << 100*(double)decayCount_phi/lineCount << ")" << std::flush;
    //
    decayCount_phi++;
  }
  std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been decayed. (" << 100.*(double)decayCount_phi/lineCount << ")" << std::endl;
  std::cout << decayCount_phi_accepted << " light dark mater particles accepted within DUNE ND fiducial volume.\n";
  std::cout << "Geometrical acceptance : " << (double)decayCount_phi_accepted / decayCount_phi * 4.0 * M_PI * dist_mc0_nd * dist_mc0_nd << "m^2" << '\n';


  for( auto& darkPhoton_index : darkPhotons_eta_dump )
  {
    decay.SetDecay(darkPhoton_index, 2, masses, "");
    decay.Generate();
    dm_particle1 = decay.GetDecay(0);
    dm_particle2 = decay.GetDecay(1);
    // TODO: determine each particle is heading toward ND fiducial volume
    //
    rayV.setX(dm_particle1->Px());
    rayV.setY(dm_particle1->Py());
    rayV.setZ(dm_particle1->Pz());
    hphi_eta_dump->Fill(dm_particle1->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle1);
      hphi_eta_dump_acc->Fill(dm_particle1->E());
      decayCount_phi_accepted++;
    }
    rayV.setX(dm_particle2->Px());
    rayV.setY(dm_particle2->Py());
    rayV.setZ(dm_particle2->Pz());
    hphi_eta_dump->Fill(dm_particle2->E());

    interCube = RayIntersectsCube(rayO, rayV, &box);
    if( interCube )
    {
      darkMatters.push_back(*dm_particle2);
      hphi_eta_dump_acc->Fill(dm_particle2->E());
      decayCount_phi_accepted++;
    }
    if( decayCount_phi % 100000 == 0 ) std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been loaded. (" << 100*(double)decayCount_phi/lineCount << ")" << std::flush;
    //
    decayCount_phi++;
  }
  std::cout << '\r' << "Decaying A's into dark matter pair ... " << decayCount_phi << " particles have been decayed. (" << 100.*(double)decayCount_phi/lineCount << ")" << std::endl;
  std::cout << decayCount_phi_accepted << " light dark mater particles accepted within DUNE ND fiducial volume.\n";
  std::cout << "Geometrical acceptance : " << (double)decayCount_phi_accepted / decayCount_phi * 4.0 * M_PI * dist_mc0_nd * dist_mc0_nd << "m^2" << '\n';

  // Calculate expected number of signals
  //


  // Producing sensitivity curve
  //

  // Scale the outputs
  //
  const double sim_pot = 1.0e+6;
  const double scaleFactor = pot_per_year / sim_pot / fEx_width;
  hpizero_target->Scale(scaleFactor);
  hpizero_dump->Scale(pot_per_year / 1.0e+4 / fEx_width);
  heta_target->Scale(scaleFactor);
  heta_dump->Scale(pot_per_year / 1.0e+4 / fEx_width);
  hAprime_pi0_target->Scale(scaleFactor);
  hAprime_pi0_dump->Scale(scaleFactor);
  hAprime_eta_target->Scale(scaleFactor);
  hAprime_eta_dump->Scale(scaleFactor);
  hphi_pi0_target->Scale(scaleFactor);
  hphi_pi0_target_acc->Scale(scaleFactor);
  hphi_pi0_dump->Scale(scaleFactor);
  hphi_pi0_dump_acc->Scale(scaleFactor);
  hphi_eta_target->Scale(scaleFactor);
  hphi_eta_target_acc->Scale(scaleFactor);
  hphi_eta_dump->Scale(scaleFactor);
  hphi_eta_dump_acc->Scale(scaleFactor);

  // Write output
  //
  TFile fOutput("proto_output.root","RECREATE");
  hpizero_target->Write();
  hpizero_dump->Write();
  heta_target->Write();
  heta_dump->Write();
  hAprime_pi0_target->Write();
  hAprime_pi0_dump->Write();
  hAprime_eta_target->Write();
  hAprime_eta_dump->Write();
  hphi_pi0_target->Write();
  hphi_pi0_target_acc->Write();
  hphi_pi0_dump->Write();
  hphi_pi0_dump_acc->Write();
  hphi_eta_target->Write();
  hphi_eta_target_acc->Write();
  hphi_eta_dump->Write();
  hphi_eta_dump_acc->Write();
  */
  /*
  TH1D* hphi_eff = (TH1D*)hphi_acc->Clone("hphi_eff");
  hphi_eff->Divide(hphi);
  hphi_eff->SetTitle("Efficiency");
  hphi_eff->GetYaxis()->SetTitle("Efficiency");
  hphi_eff->Write();
  */

  //fOutput.Close();

	return 0;
}

