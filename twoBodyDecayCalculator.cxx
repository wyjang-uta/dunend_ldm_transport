#include "twoBodyDecayCalculator.h"

/// @brief default constructor
twoBodyDecayCalculator::twoBodyDecayCalculator()
 : f_mother_mass(-1),
   decay_counter(0)
{
    randomSeed = 12345;
}

/// @brief constructor
/// @param mass             mass of mother particle
/// @param mass_daughter1   mass of first daughter particle
/// @param mass_daughter2   mass of second daughter particle
twoBodyDecayCalculator::twoBodyDecayCalculator(double mass, double mass_daughter1, double mass_daughter2) 
  : f_mother_mass(mass),
    f_daughter_masses{mass_daughter1, mass_daughter2},
    decay_counter(0)

{
    randomSeed = 12345;
}

/// @brief constructor with file name specified
/// @param[in] name of file to load.
twoBodyDecayCalculator::twoBodyDecayCalculator(const std::string& inputFileName, double mass, double mass_daughter1, double mass_daughter2)
 : f_mother_mass(mass),
   f_daughter_masses{mass_daughter1, mass_daughter2},
   decay_counter(0)
{
    randomSeed = 12345;
    f_fileName = inputFileName;
    f_fileStream.open(f_fileName);
    const size_t bufferSize = 1024 * 1024; // size of the buffer (1 MB)
    std::vector<char> buffer(1024*1024);
    f_fileStream.rdbuf()->pubsetbuf(buffer.data(),bufferSize);
}

twoBodyDecayCalculator::twoBodyDecayCalculator(const std::vector<TLorentzVector>& mother, double mass, double mass_daughter1, double mass_daughter2)
 : f_mother_mass(mass),
   f_daughter_masses{mass_daughter1, mass_daughter2},
   decay_counter(0)
{
   f_mothers = mother;
}

/// @brief destructor
twoBodyDecayCalculator::~twoBodyDecayCalculator() {
    if( f_fileStream.is_open() ) {
        f_fileStream.close();
    }
}

/// @brief return particles vector container
/// @return this method returns the vector container holding the 4-vector information of mother particles that to be decayed.
const std::vector<TLorentzVector>& twoBodyDecayCalculator::getMotherParticles() const {
    return f_mothers;
}


/// @brief return daughter particles vector container
/// @return this method returns the vector container of daughter particles. the
///         function check whether the decay operation had happened or not and 
///         when decay operation happened, then return the vector container 
///         that holds 4-vectors of daughter particles
const std::vector<std::pair<int, TLorentzVector>>& twoBodyDecayCalculator::getDaughterParticles() const {
    return f_daughters;
}

/// @brief  return the vector containter of the first decay product
/// @return 
const std::vector<TLorentzVector>& twoBodyDecayCalculator::getDaughter1() const {
    return f_daughter1;
}

/// @brief  return the vector containter of the second decay product
/// @return 
const std::vector<TLorentzVector>& twoBodyDecayCalculator::getDaughter2() const {
    return f_daughter2;
}

/// @brief return the filestream instance
/// @return filestream instance
const std::ifstream& twoBodyDecayCalculator::getFileStream() const {
    return f_fileStream;
}

/// @brief set the filestream
/// @param f_fileName 
/// @return file status boolean flag
bool twoBodyDecayCalculator::setFileStream(const std::string& f_fileName) {
    f_fileStream.open(f_fileName);
    return f_fileStream.is_open();
}

const double twoBodyDecayCalculator::getMass() const {
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
    int progress_interval = 1000000;
    unsigned long long currentProgress = 0;
    unsigned long long newProgress = 0;

    TLorentzVector particle;
    while( f_fileStream >> px >> py >> pz >> E ) {
        particle.SetPxPyPzE(px, py, pz, E);
        f_mothers.push_back(particle);

        // progress indicator
        newProgress = f_mothers.size();
        if( newProgress % progress_interval == 0 )
        {
            currentProgress = newProgress;
            std::cout << '\r' << "Loading data from " << f_fileName << " ... [" << currentProgress << "] ... " << progress_symbols[progress_symbol_index] << std::flush;
            progress_symbol_index = ( progress_symbol_index + 1 ) % 4;
        }
    }

    std::cout << '\r' << "Loading data from " << f_fileName << " ... [DONE] -- " << f_mothers.size() << " particles loaded\n";

    return f_mothers.size();
}

int twoBodyDecayCalculator::unload() {
    f_mothers.clear();
    return 0;
}

unsigned int twoBodyDecayCalculator::decay() {
    unsigned int data_size = f_mothers.size();
    if( data_size == 0 ) {
        std::cerr<< "no data have been loaded.\n";
        return f_mothers.size();
    }
    char progress_symbols[4] = {'-', '/', '|', '\\'};
    int progress_symbol_index = 0;
    int progress_interval = 1000000;
    unsigned long long currentProgress = 0;
    unsigned long long newProgress = 0;

    unsigned int count = 0;
    for( auto& mother : f_mothers ) {
        f_decay.SetDecay(mother, 2, f_daughter_masses, "");
        f_decay.Generate();
        f_daughter1.push_back(*f_decay.GetDecay(0));
        f_daughter2.push_back(*f_decay.GetDecay(1));
        f_daughters.push_back(std::make_pair(0, *(f_decay.GetDecay(0))));
        f_daughters.push_back(std::make_pair(1, *(f_decay.GetDecay(1))));

        newProgress = count;
        if( newProgress % progress_interval == 0 )
        {
            currentProgress = newProgress;
            std::cout << '\r' << "Calculating 2-body decay ... [" << currentProgress << "] ... " << progress_symbols[progress_symbol_index] << std::flush;
            progress_symbol_index = ( progress_symbol_index + 1 ) % 4;
        }

        count++;
    }
    std::cout << '\r' << "Calculating 2-body decay ... [ DONE ]                          \n";

    decay_counter = count;

    return count;
}

const unsigned int twoBodyDecayCalculator::getDecayCount() const {
    return decay_counter;
}
