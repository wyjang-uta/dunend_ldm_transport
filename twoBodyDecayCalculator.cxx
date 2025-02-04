#include "twoBodyDecayCalculator.h"

namespace RTMath = ROOT::Math;

/// @brief default constructor
twoBodyDecayCalculator::twoBodyDecayCalculator() : particleCounter(0) {}

/// @brief constructor with file name specified
/// @param[in] name of file to load.
twoBodyDecayCalculator::twoBodyDecayCalculator(const std::string& inputFileName) : particleCounter(0) {
    fileName = inputFileName;
    fileStream.open(fileName);
}

/// @brief constructor with file name specified
/// @param[in] name of file to load.
twoBodyDecayCalculator::twoBodyDecayCalculator(const char* inputFileName) : particleCounter(0) {
    fileName.append(inputFileName);
    fileStream.open(fileName);
}

/// @brief destructor
twoBodyDecayCalculator::~twoBodyDecayCalculator() {
    if( fileStream.is_open() ) {
        fileStream.close();
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
    return fileStream;
}

/// @brief set the filestream
/// @param fileName 
/// @return file status boolean flag
bool twoBodyDecayCalculator::setFileStream(const std::string& fileName) {
    fileStream.open(fileName);
    return fileStream.is_open();
}

/// @brief load the data from the specified file.
/// @param mothers 
/// @return 
int twoBodyDecayCalculator::load() {
    if( !fileStream.is_open() ) {
        std::cerr << "[twoBodyDecayCalculator]: File (" << fileName << ")was not opened." << std::endl;
        exit(EXIT_FAILURE);
    }

    double px, py, pz, E;
    char progress_symbols[4] = {'-', '/', '|', '\\'};
    int progress_symbol_index = 0;
    unsigned long long currentProgress = 0;
    unsigned long long newProgress = 0;

    RTMath::LorentzVector<RTMath::PxPyPzE4D<double>> particle;
    while( fileStream >> px >> py >> pz >> E ) {
        particle.SetPxPyPzE(px, py, pz, E);
        mothers.push_back(particle);
        particleCounter++;

        // progress indicator
        newProgress = particleCounter;
        if( newProgress > currentProgress )
        {
            currentProgress = newProgress;
            std::cout << '\r' << "Loading data from " << fileName << " ... [" << currentProgress << "] ... " << progress_symbols[progress_symbol_index] << std::flush;
            progress_symbol_index = ( progress_symbol_index + 1 ) % 4;
        }
    }

    return currentProgress;
}

int twoBodyDecayCalculator::unload() {
    mothers.clear();
    return 0;
}
