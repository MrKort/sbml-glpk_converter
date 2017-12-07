#include <chrono>
#include "bacterium.h"

int main (int argc, char* argv[])
{
    try {
        //Start of sim timer
        auto startTime = std::chrono::high_resolution_clock::now();

        if(!(argc > 1)){
            throw std::string ("no arguments (argc/argv) were given...");
        }

        int bacNum = argc - 1;
        std::vector<const char*> inputFiles(bacNum);
        for(int i = 0; i < bacNum; ++i) {
            inputFiles[i] = argv[i + 1];
        }

        std::vector<Bacterium> bac(bacNum);

        for(int i = 0; i < bacNum; ++i) {
            bac[i].readInputSBML(inputFiles[i]);
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
