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

namespace RTMath = ROOT::Math;

class twoBodyDecayCalculator
{
    public:
        twoBodyDecayCalculator();
        twoBodyDecayCalculator(const std::string& fileName);
        twoBodyDecayCalculator(const char* fileName);
        ~twoBodyDecayCalculator();

        std::vector<RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>& getMotherParticles();
        std::vector<std::pair<int, RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>>& getDaughterParticles();
        std::ifstream& getFileStream();
        bool setFileStream(const std::string& fileName);

        int load();
        int unload();

        int particleCounter;

    private:
        std::string fileName;
        std::ifstream fileStream;

        std::vector<RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>> mothers;
        std::vector<std::pair<int, RTMath::LorentzVector<RTMath::PxPyPzE4D<double>>>> daughters;
};

#endif
