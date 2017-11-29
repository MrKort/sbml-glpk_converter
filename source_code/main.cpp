#include <chrono>

#include "bacterium.h"

int main (int argc, char* argv[])
{
    //Start of sim timer
    auto startTime = std::chrono::high_resolution_clock::now();

    Bacterium test; //Use constructor with argv to initialise
    Bacterium test2;

    //Declare test input file
    std::vector<const char*> inputFiles(2);
    inputFiles[0]= "../test.xml";
    inputFiles[1]= "../test2.xml";

    test.readInputSBML(inputFiles[0]);
    test.doFBA();
    test.outputMatrix();

    std::cout << "###################################################################\n\n";

    test2.readInputSBML(inputFiles[1]);
    test2.doFBA();
    test2.outputMatrix();

    //End of sim timer
    auto endTime = std::chrono::high_resolution_clock::now();
    double simulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 1000.0;
    std::cout << "Simulation has finished in " << simulationTime << " seconds.\n";

    return 0;
}
