#include "twoBodyDecayCalculator.h"

namespace RTMath = ROOT::Math;

/// @brief default constructor
twoBodyDecayCalculator::twoBodyDecayCalculator()
 : particleCounter(0), 
   f_mother_mass(-1)
{
}

twoBodyDecayCalculator::twoBodyDecayCalculator(double mass, double mass_daughter1, double mass_daughter2) 
  : particleCounter(0), 
    f_mother_mass(mass),
    f_daughter_masses{mass_daughter1, mass_daughter2}
{
}

/// @brief constructor with file name specified
/// @param[in] name of file to load.
twoBodyDecayCalculator::twoBodyDecayCalculator(const std::string& inputFileName, double mass, double mass_daughter1, double mass_daughter2)
 : particleCounter(0),
   f_mother_mass(mass),
   f_daughter_masses{mass_daughter1, mass_daughter2}
{
    f_fileName = inputFileName;
    f_fileStream.open(f_fileName);
}

/// @brief destructor
twoBodyDecayCalculator::~twoBodyDecayCalculator() {
    if( f_fileStream.is_open() ) {
        f_fileStream.close();
    }
}

/// @brief return particles vector container
/// @return this method returns the vector container holding the 4-vector information of mother particles that to be decayed.
std::vector<RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>& twoBodyDecayCalculator::getMotherParticles() {
    return mothers;
}


/// @brief return daughter particles vector container
/// @return this method returns the vector container of daughter particles. the
///         function check whether the decay operation had happened or not and 
///         when decay operation happened, then return the vector container 
///         that holds 4-vectors of daughter particles
std::vector<std::pair<int, RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>>& twoBodyDecayCalculator::getDaughterParticles() {
    return daughters;
}

/// @brief return the filestream instance
/// @return filestream instance
std::ifstream& twoBodyDecayCalculator::getFileStream() {
    return f_fileStream;
}

/// @brief set the filestream
/// @param f_fileName 
/// @return file status boolean flag
bool twoBodyDecayCalculator::setFileStream(const std::string& f_fileName) {
    f_fileStream.open(f_fileName);
    return f_fileStream.is_open();
}

double twoBodyDecayCalculator::getMass() {
    return f_mother_mass;
}

/// @brief load the data from the specified file.
/// @param mothers 
/// @return 
int twoBodyDecayCalculator::load() {
    if( !f_fileStream.is_open() ) {
        std::cerr << "[twoBodyDecayCalculator]: File (" << f_fileName << ")was not opened." << std::endl;
        exit(EXIT_FAILURE);
    }

    double px, py, pz, E;
    char progress_symbols[4] = {'-', '/', '|', '\\'};
    int progress_symbol_index = 0;
    unsigned long long currentProgress = 0;
    unsigned long long newProgress = 0;

    RTMath::LorentzVector<RTMath::PxPyPzE4D<double>> particle;
    while( f_fileStream >> px >> py >> pz >> E ) {
        particle.SetPxPyPzE(px, py, pz, E);
        mothers.push_back(particle);
        particleCounter++;

        // progress indicator
        newProgress = particleCounter;
        if( newProgress > currentProgress )
        {
            currentProgress = newProgress;
            std::cout << '\r' << "Loading data from " << f_fileName << " ... [" << currentProgress << "] ... " << progress_symbols[progress_symbol_index] << std::flush;
            progress_symbol_index = ( progress_symbol_index + 1 ) % 4;
        }
    }

    return currentProgress;
}

int twoBodyDecayCalculator::unload() {
    mothers.clear();
    return 0;
}

int twoBodyDecayCalculator::decay() {
    unsigned long long decayCounter = 0;
    if( particleCounter == 0 ) {
        std::cerr<< "no data have been loaded.\n";
        return decayCounter;
    }

    return decayCounter;
}