#ifndef BACTERIUM_H
#define BACTERIUM_H

#include <glpk.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/fbcfwd.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
#include <map>
#include <fstream>

LIBSBML_CPP_NAMESPACE_USE

class Bacterium {
public:
    void readInputSBML(const char *&);
    void doFBA();
    void outputMatrix();

private:
    //Species ID from SBML file
    std::string speciesID;

    //Index associated with the objective function reaction (biomass production)
    unsigned int ObjFunc;
    unsigned int revCount;

    //Declare internal and external maps to link SBML metabolite name to GLPK metabolite index
    std::map<std::string, int> environment;
    std::map<std::string, int> cytosol;

    //Declare metabolite and reaction maps to link GLPK metabolite index to SBML metabolite name
    std::map<int, std::string> metaName;
    std::map<int, std::string> reacName;

    //Delare vectors to be used by GLPK
    std::vector<int> vecMeta, vecReac;          //vector with GLPK row and column index values
    std::vector<double> vecStoi;                //vector with GLPK stoichiometry values
};


void Bacterium::readInputSBML(const char *&inputFile){
//##############################//
//  SBML FILE READING FUNCTION  //
//##############################//

        //Read test input file with SBML
        SBMLDocument *document = readSBML(inputFile);

        unsigned int errors = document->getNumErrors();
        if(errors){
            std::string errorMessage = std::string("could not read SBML file: ") + inputFile;
            throw std::string(errorMessage);
        }

        Model *model = document->getModel();
        speciesID = model->getId();

        ListOfSpecies *spList = model->getListOfSpecies();
        ListOfReactions *reacList = model->getListOfReactions();

        //Read metabolite data from SBML file
        for(unsigned int i = 0u, metaCount = 0u; i < spList->size(); ++i) {
            std::string spId = spList->get(i)->getId();
            std::string spComp = spId.substr(spId.size()-1, spId.size());

            //Store metabolite name and number in map
            if(spComp == "e") {
                metaName[metaCount] = spId; //compartmentalisation independent storage map
                environment[spId] = metaCount++;
            }
            else if(spComp == "c") {
                metaName[metaCount] = spId; //compartmentalisation independent storage map
                cytosol[spId] = metaCount++;
            }
            else {
                std::string errorMessage = std::string("unexpected compartment whilst reading SBML file: ") + inputFile;
                throw std::string (errorMessage);
            }
        }
        if(!cytosol.find("M_biomass_c")->second) {
            throw std::string ("could not find a biomass metabolite in cytosol compartment.");
        }

        //These counters are initialised at 1, because glp can't handle 0's
        unsigned int matrixCount = 1u;
        revCount = 1u;

        //Write substrates and products from SMBL files into GLPK arrays
        for(unsigned int i = 0u; i < reacList->size(); ++i) {
            const Reaction *reac = reacList->get(i);
            const std::string reacId = reac->getId();
            reacName[i] = reacId;

            //Determine wheter the reaction name is biomass0, this is usually the objective function
            //In the future update the objective function determenation by using SBML::FBC:objective_function
            if(reac->getName() == "biomass0") {
                ObjFunc = i;
            }

            //For loop for substrates
            unsigned int m = reac->getListOfReactants()->size();
            for(unsigned int j = 0u; j < m; ++j) {
                const std::string subId = reac->getReactant(j)->getSpecies();
                const std::string subComp = subId.substr(subId.size()-1, subId.size());
                const double subStoi = reac->getReactant(j)->getStoichiometry();

                //Build stoichiometry matrix using maps.find() to get metabolite number
                if(subComp == "e") {
                    const int subNum = environment.find(subId)->second;
                    vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                    ++matrixCount;
                    if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                        vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
                        ++matrixCount;
                    }
                }
                else if(subComp == "c") {
                    const int subNum = cytosol.find(subId)->second;
                    vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                    ++matrixCount;
                    if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                        vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
                        ++matrixCount;
                    }
                }
                else {
                    std::string errorMessage = std::string("unexpected compartment whilst reading substrates from SBML file: ") + inputFile;
                    throw std::string (errorMessage);
                }
            }

            //For loop for products
            m = reac->getListOfProducts()->size();
            for(unsigned int j = 0u; j < m; ++j) {
                const std::string prodId = reac->getProduct(j)->getSpecies();
                const std::string prodComp = prodId.substr(prodId.size()-1, prodId.size());
                const double prodStoi = reac->getProduct(j)->getStoichiometry();

                //Build stoichiometry matrix using maps.find() to get metabolite number
                if(prodComp == "e") {
                    const int prodNum = environment.find(prodId)->second;
                    vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                    ++matrixCount;
                    if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                        vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                        ++matrixCount;
                    }
                }
                else if(prodComp == "c") {
                    const int prodNum = cytosol.find(prodId)->second;
                    vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                    ++matrixCount;
                    if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                        vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                        ++matrixCount;
                    }
                }
                else {
                    std::string errorMessage = std::string("unexpected compartment whilst reading products from SBML file: ") + inputFile;
                    throw std::string (errorMessage);
                }
            }
            if(reac->getReversible()) {
                reacName[reacList->size() + revCount - 1] = reacId + "_reverse";
                ++revCount;
            }
        }

        //Done reading file, now delete the SBML document variable
        delete document;
        std::cout << "Succesfully read data from SBML input file: " << inputFile
                  << "\nMetabolic data for species: " << speciesID << " stored.\n\n";
}

void Bacterium::doFBA(){
//##############################//
// GLPK READ AND SOLVE FUNCTION //
//##############################//

//CREATE A NEW GLPK PROBLEM
            glp_prob *lp;
            lp = glp_create_prob();
            glp_set_obj_dir(lp, GLP_MAX);

//ADD ROW AND COLUMNS FOR STOICHIOMETRY MATRIX
            glp_add_cols(lp, reacName.size());
            glp_add_rows(lp, metaName.size());

//ASSIGN METABOLITE NAME AND BOUNDS
            for(unsigned int i = 0; i < metaName.size(); ++i){
                glp_set_row_name(lp, i+1, metaName.find(i)->second.c_str());
                glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
            }

//ASSIGN REACTION NAME AND BOUNDS

            //Need to update this depending on the environment...
            //Separate from first GLP initialisation

            for(unsigned int i = 0; i < reacName.size(); ++i) {
                glp_set_col_name(lp, i+1, reacName.find(i)->second.c_str());
                glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, 1000.0);

                if(i == ObjFunc) { //check for objective function reaction
                    glp_set_obj_coef(lp, i+1, 1.0);
                }
                else {
                    glp_set_obj_coef(lp, i+1, 0.0);
                }
            }

//LOAD AND SOLVE GLPK PROBLEM
            //vector.data() to access the internal arrays
            //vector.data()-1 to skip the first[0] element
            glp_load_matrix(lp, vecStoi.size(),vecMeta.data()-1,vecReac.data()-1,vecStoi.data()-1);
            glp_simplex(lp, NULL);

//DELETE THE GLPK PROBLEM

            //Is deleting the problem useful?
            //Only need to update column bounds depending on environment...

            glp_delete_prob(lp);
            std::cout << "Succesfully finished Linear Programming step.\n\n";
}

void Bacterium::outputMatrix(){
//##############################//
//    TEST OUTPUT FUNCTIONS     //
//##############################//

        std::cout << "Number of external metabolites: " << environment.size() << "\n";
        std::cout << "Number of internal metabolites: " << cytosol.size() << "\n";
        std::cout << "Total number metabolites: " << metaName.size() << "\n";
        std::cout << "\nNumber of forward reactions: " << reacName.size() - (revCount-1) << "\n";
        std::cout << "Number of reverse reactions: " << revCount - 1 << "\n";
        std::cout << "Total number of reactions: " << reacName.size() << "\n";
        std::cout << "\nNumber of non-zeroes in sparse matrix: " << vecStoi.size() << "\n\n";

        //Output sparse matrix with reaction and metabolite IDs to .csv file
        std::string outputFile = "./Sparse_Matrices/SM_" + speciesID + ".csv";
//Check if Sparse_Matrices/ folder exists, otherwise make it
        std::ofstream matrix(outputFile);
        if(!matrix.is_open()) {
            std::string errorMessage = std::string("could not open outputfile: ") + outputFile;
            throw std::string (errorMessage);
        }

        matrix << "Matrix Index\tMetabolite#\tMetaboliteID\tReaction#\tReactionID\tStoichiometry\n";
        for(unsigned int i=1; i < vecStoi.size(); ++i){
            matrix << i << "\t" << vecMeta[i-1] << "\t" << metaName.find(vecMeta[i-1]-1)->second
                        << "\t" << vecReac[i-1] << "\t" << reacName.find(vecReac[i-1]-1)->second
                        << "\t" << vecStoi[i-1] << "\n";
        }
        matrix.close();
        std::cout << "Succesfully created sparse matrix in file: " << outputFile << "\n\n";
}

#endif // BACTERIUM_H
