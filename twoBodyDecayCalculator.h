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
        twoBodyDecayCalculator(const std::vector<TLorentzVector>& mothers, double mass, double mass_daughter1, double mass_daughter2);
        ~twoBodyDecayCalculator();

        const std::vector<TLorentzVector>& getMotherParticles() const;
        const std::vector<std::pair<int, TLorentzVector>>& getDaughterParticles() const;
        const std::vector<TLorentzVector>& getDaughter1() const;
        const std::vector<TLorentzVector>& getDaughter2() const;
        const std::ifstream& getFileStream() const;
        bool setFileStream(const std::string& fileName);
        const double getMass() const;
        const unsigned int getDecayCount() const;

        // 주요 함수들
        int load();
        int unload();
        unsigned int decay();


    private:
        std::string f_fileName;
        std::ifstream f_fileStream;

        double f_mother_mass;
        double f_daughter_masses[2];

        TGenPhaseSpace f_decay;
        std::vector<TLorentzVector> f_mothers;
        std::vector<TLorentzVector> f_daughter1;
        std::vector<TLorentzVector> f_daughter2;
        std::vector<std::pair<int, TLorentzVector>> f_daughters;

        int randomSeed;
        unsigned int decay_counter;
};

#endif
