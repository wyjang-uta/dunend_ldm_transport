#include "twoBodyDecayCalculator.h"

/// @brief default constructor
twoBodyDecayCalculator::twoBodyDecayCalculator()
 : f_mother_mass(-1)
{
}

twoBodyDecayCalculator::twoBodyDecayCalculator(double mass, double mass_daughter1, double mass_daughter2) 
  : f_mother_mass(mass),
    f_daughter_masses{mass_daughter1, mass_daughter2}
{
}

/// @brief constructor with file name specified
/// @param[in] name of file to load.
twoBodyDecayCalculator::twoBodyDecayCalculator(const std::string& inputFileName, double mass, double mass_daughter1, double mass_daughter2)
 : f_mother_mass(mass),
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
std::vector<TLorentzVector>& twoBodyDecayCalculator::getMotherParticles() {
    return f_mothers;
}


/// @brief return daughter particles vector container
/// @return this method returns the vector container of daughter particles. the
///         function check whether the decay operation had happened or not and 
///         when decay operation happened, then return the vector container 
///         that holds 4-vectors of daughter particles
std::vector<std::pair<int, TLorentzVector>>& twoBodyDecayCalculator::getDaughterParticles() {
    return f_daughters;
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
        std::cerr << "[twoBodyDecayCalculator]: File (" << f_fileName << ") was not opened." << std::endl;
        exit(EXIT_FAILURE);
    }

    double px, py, pz, E;
    char progress_symbols[4] = {'-', '/', '|', '\\'};
    int progress_symbol_index = 0;
    unsigned long long currentProgress = 0;
    unsigned long long newProgress = 0;

    TLorentzVector particle;
    while( f_fileStream >> px >> py >> pz >> E ) {
        particle.SetPxPyPzE(px, py, pz, E);
        f_mothers.push_back(particle);

        // progress indicator
        newProgress = f_mothers.size();
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
    f_mothers.clear();
    return 0;
}

int twoBodyDecayCalculator::decay() {
    unsigned long long decayCounter = 0;
    if( f_mothers.size() == 0 ) {
        std::cerr<< "no data have been loaded.\n";
        return decayCounter;
    }

    for( auto& mother : f_mothers ) {
        f_decay.SetDecay(mother, 2, f_daughter_masses, "");
        f_decay.Generate();
    }


    return decayCounter;
}
