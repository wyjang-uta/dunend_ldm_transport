#ifndef TWOBODYDECAYCALCULATOR_H
#define TWOBODYDECAYCALCULATOR_H

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <utility>

#include "Math/LorentzVector.h"
#include "TGenPhaseSpace.h"

namespace RTMath = ROOT::Math;

class twoBodyDecayCalculator
{
    public:
        twoBodyDecayCalculator();
        twoBodyDecayCalculator(double mass, double mass_daughter1, double mass_daughter2);
        twoBodyDecayCalculator(const std::string& fileName, double mass, double mass_daughter1, double mass_daughter2);
        ~twoBodyDecayCalculator();

        std::vector<RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>& getMotherParticles();
        std::vector<std::pair<int, RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>>& getDaughterParticles();
        std::ifstream& getFileStream();
        bool setFileStream(const std::string& fileName);
        double getMass();

        int load();
        int unload();
        int decay();

        int particleCounter;

    private:
        std::string f_fileName;
        std::ifstream f_fileStream;

        double f_mother_mass;
        double f_daughter_masses[2];

        TGenPhaseSpace f_decay;
        std::vector<RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>> mothers;
        std::vector<std::pair<int, RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>> daughters;
};

#endif
