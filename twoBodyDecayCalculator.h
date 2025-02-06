#ifndef TWOBODYDECAYCALCULATOR_H
#define TWOBODYDECAYCALCULATOR_H

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <utility>

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

class twoBodyDecayCalculator
{
    public:
        twoBodyDecayCalculator();
        twoBodyDecayCalculator(double mass, double mass_daughter1, double mass_daughter2);
        twoBodyDecayCalculator(const std::string& fileName, double mass, double mass_daughter1, double mass_daughter2);
        ~twoBodyDecayCalculator();

        std::vector<TLorentzVector>& getMotherParticles();
        std::vector<std::pair<int, TLorentzVector>>& getDaughterParticles();
        std::ifstream& getFileStream();
        bool setFileStream(const std::string& fileName);
        double getMass();

        // 주요 함수들
        int load();
        int unload();
        int decay();


    private:
        std::string f_fileName;
        std::ifstream f_fileStream;

        double f_mother_mass;
        double f_daughter_masses[2];

        TGenPhaseSpace f_decay;
        std::vector<TLorentzVector> f_mothers;
        std::vector<std::pair<int, TLorentzVector>> f_daughters;
};

#endif
