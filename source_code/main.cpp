#include <chrono>
#include "bacterium.h"

int main (int argc, char* argv[])
{
    try {
//Start of sim timer
        auto startTime = std::chrono::high_resolution_clock::now();

//Check if arguments are given
        if(!(argc > 1)){
            throw std::string ("please provide at least one SBML file as argument...");
        }

//Put arguments into inputFiles vector
        int bacNum = argc - 1;
        std::vector<const char*> inputFiles(bacNum);
        for(int i = 0; i < bacNum; ++i) {
            inputFiles[i] = argv[i + 1];
        }

//Initialise a vector of Bacteria (size determined by the amount of arguments)
        std::vector<Bacterium> bac(bacNum);

//Read SBML and do GLPK for each Bacteria
        for(int i = 0; i < bacNum; ++i) {
            bac[i].readFileSBML(inputFiles[i]);
            bac[i].doFBA();
            bac[i].outputMatrix();
            if (i != bacNum - 1) {
                std::cout << "###################################################################\n\n";
            }
        }

//End of sim timer
        auto endTime = std::chrono::high_resolution_clock::now();
        double simulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 1000.0;
        std::cout << "Simulation has finished in " << simulationTime << " seconds.\n";
    }

    catch(std::string &errorMessage) {
        std::cerr << "Encountered an error: " << errorMessage << '\n';
    }
    return 0;
}
