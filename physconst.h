// all constants refered from the Review of Particle Physics 2024
// all constants in natural unit system
const double alpha                    = 0.0072973525693;                  // 7.297 352 5693(11)e-3 = 1/137.035 999 084(21), uncertainty is 0.15 in ppb(1e-9).
const double alpha_D                  = 0.1;                              // dark fine structure constant
const double hbarc                    = 197.3269804;                      // 197.326 980 4 MeV*fm
const double epsilon                  = 0.01;                             // kinetic mixing parameter
const double four_pi_epsilon_squared  = 0.0012566371;                     // 4*pi*epsilon^2 in natural unit when epsilon = 0.01
const double avogadro                 = 6.022e+23;

// particle masses
const double m_e                      = 0.51099895000;                    // 0.510 998 950 00(15) MeV/c2 = 9.109 383 7015(28)e-31 kg, with uncertainty 0.30 ppb(1e-9)
const double m_pi0                    = 134.9768;                         // 134.9768 +/- 0.0005 MeV/c2

// branching ratios
const double br_eta_gammagamma        = 0.3936;                           // 39.36 +/- 0.18 %

// DUNE near detector dimensions (lengths in m)
const double active_x                 = 7.0;
const double active_y                 = 3.0;
const double active_z                 = 5.0;

const double fiducial_x               = 6.0;
const double fiducial_y               = 2.0;
const double fiducial_z               = 3.0;

const double fiducial_vertex[8][3] = {
  { 0.5 * fiducial_x,  0.5 * fiducial_y, -0.5 * fiducial_z},
  {-0.5 * fiducial_x,  0.5 * fiducial_y, -0.5 * fiducial_z},
  { 0.5 * fiducial_x, -0.5 * fiducial_y, -0.5 * fiducial_z},
  {-0.5 * fiducial_x, -0.5 * fiducial_y, -0.5 * fiducial_z},
  { 0.5 * fiducial_x,  0.5 * fiducial_y,  0.5 * fiducial_z},
  {-0.5 * fiducial_x,  0.5 * fiducial_y,  0.5 * fiducial_z},
  { 0.5 * fiducial_x, -0.5 * fiducial_y,  0.5 * fiducial_z},
  {-0.5 * fiducial_x, -0.5 * fiducial_y,  0.5 * fiducial_z}
};

const double dist_mc0_nd              = 574.0;
const double dist_ground_nd           = 62.50;
const double dune_nd_rotation_angle_in_x = asin(dist_ground_nd/dist_mc0_nd);

const double cosine = cos(dune_nd_rotation_angle_in_x);
const double sine   = sin(dune_nd_rotation_angle_in_x);
const double rotated_fiducial_vertex[8][3] = {
  { fiducial_vertex[0][0], cosine * fiducial_vertex[0][1] - sine * fiducial_vertex[0][2], sine * fiducial_vertex[0][1] + cosine * fiducial_vertex[0][2]},
  { fiducial_vertex[1][0], cosine * fiducial_vertex[1][1] - sine * fiducial_vertex[1][2], sine * fiducial_vertex[1][1] + cosine * fiducial_vertex[1][2]},
  { fiducial_vertex[2][0], cosine * fiducial_vertex[2][1] - sine * fiducial_vertex[2][2], sine * fiducial_vertex[2][1] + cosine * fiducial_vertex[2][2]},
  { fiducial_vertex[3][0], cosine * fiducial_vertex[3][1] - sine * fiducial_vertex[3][2], sine * fiducial_vertex[3][1] + cosine * fiducial_vertex[3][2]},
  { fiducial_vertex[4][0], cosine * fiducial_vertex[4][1] - sine * fiducial_vertex[4][2], sine * fiducial_vertex[4][1] + cosine * fiducial_vertex[4][2]},
  { fiducial_vertex[5][0], cosine * fiducial_vertex[5][1] - sine * fiducial_vertex[5][2], sine * fiducial_vertex[5][1] + cosine * fiducial_vertex[5][2]},
  { fiducial_vertex[6][0], cosine * fiducial_vertex[6][1] - sine * fiducial_vertex[6][2], sine * fiducial_vertex[6][1] + cosine * fiducial_vertex[6][2]},
  { fiducial_vertex[7][0], cosine * fiducial_vertex[7][1] - sine * fiducial_vertex[7][2], sine * fiducial_vertex[7][1] + cosine * fiducial_vertex[7][2]}
};

const double nd_location_vertex[8][3] = {
  { rotated_fiducial_vertex[0][0], rotated_fiducial_vertex[0][1], rotated_fiducial_vertex[0][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[1][0], rotated_fiducial_vertex[1][1], rotated_fiducial_vertex[1][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[2][0], rotated_fiducial_vertex[2][1], rotated_fiducial_vertex[2][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[3][0], rotated_fiducial_vertex[3][1], rotated_fiducial_vertex[3][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[4][0], rotated_fiducial_vertex[4][1], rotated_fiducial_vertex[4][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[5][0], rotated_fiducial_vertex[5][1], rotated_fiducial_vertex[5][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[6][0], rotated_fiducial_vertex[6][1], rotated_fiducial_vertex[6][2] + dist_mc0_nd },
  { rotated_fiducial_vertex[7][0], rotated_fiducial_vertex[7][1], rotated_fiducial_vertex[7][2] + dist_mc0_nd }
};

