#ifndef BACTERIUM_H
#define BACTERIUM_H

#include <glpk.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/fbcfwd.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
#include <map>
#include <fstream>
#include <exception>

LIBSBML_CPP_NAMESPACE_USE

class Bacterium {
//##############################//
//      CLASS DEFINITION        //
//##############################//
public:
    void readFileSBML(const char *&);
    void doFBA();
    void outputMatrix();

private:
    //Species ID from SBML file
    std::string speciesID;

    //Index associated with the objective function reaction (biomass production)
    unsigned int ObjFunc;

    //Total amount of forward and reverse reactions
    unsigned int fwdCount;
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

    //void functions used by readFileSBML()
    void readMetaSBML(const ListOfSpecies *&);
    void readReacSBML(const ListOfReactions *&);
    void initMatrix(const unsigned int &, const bool &, const std::string &, const std::string &, const double &);
};

//#################################################################################################

void Bacterium::readFileSBML(const char *&inputFile) {
//##############################//
//  SBML FILE READING FUNCTION  //
//##############################//

//Open SBML input file
        SBMLDocument *document = readSBML(inputFile);
//Check for errors in SBML file
        unsigned int errors = document->getNumErrors();
        if(errors){
            std::string errorMessage = std::string("could not read SBML file: ") + inputFile;
            throw std::runtime_error(errorMessage.c_str());
        }

//Extract SBML model from file
        Model *model = document->getModel();
        speciesID = model->getId();

//Read metabolite data from SBML model
        const ListOfSpecies *spList = model->getListOfSpecies();
        readMetaSBML(spList);

//Check if biomass is present as a cytosolic metabolite
        if(!cytosol.find("M_biomass_c")->second) {
            throw std::logic_error("could not find a biomass metabolite in cytosol compartment.");
        }

//Read reaction data from SBML model & generate initial GLPK matrix
        const ListOfReactions *reacList = model->getListOfReactions();
        fwdCount = reacList->size();
        revCount = 1u;
        readReacSBML(reacList);

//Delete SBML document variable
        delete document;
        std::cout << "Succesfully read data from SBML input file: " << inputFile
                  << "\nMetabolic data for species: " << speciesID << " stored.\n\n";
}

void Bacterium::readMetaSBML(const ListOfSpecies *&spList) {
//###############################//
// SBML METABOLITE READ FUNCTION //
//###############################//
    for(unsigned int metaCount = 0u; metaCount < spList->size(); ++metaCount) {
        std::string spId = spList->get(metaCount)->getId();
        std::string spComp = spId.substr(spId.size()-1, spId.size());

//Store metabolite name and number in environment or cytosol map
        if(spComp == "e") {
            metaName[metaCount] = spId; //compartmentalisation independent storage map
            environment[spId] = metaCount;
        }
        else if(spComp == "c") {
            metaName[metaCount] = spId; //compartmentalisation independent storage map
            cytosol[spId] = metaCount;
        }
        else {
            std::string errorMessage = std::string("unexpected compartment whilst reading SMBL file for species: ") + speciesID;
            throw std::logic_error(errorMessage);
        }
    }
}

void Bacterium::readReacSBML(const ListOfReactions *&reacList) {
//##############################//
// SBML REACTION READ FUNCTION  //
//##############################//
    for(unsigned int reacCount = 0u; reacCount < fwdCount; ++reacCount) {
        const Reaction *reac = reacList->get(reacCount);
        const bool reacRev = reac->getReversible();
        const std::string reacId = reac->getId();
        reacName[reacCount] = reacId;

//Determine whether the reaction name is biomass0, this is usually the objective function
        //In the future update the objective function determination by using SBML::FBC:objective_function
        if(reac->getName() == "biomass0") {
            ObjFunc = reacCount; //store index number of the "objective function"
        }

//For loop for substrates
        unsigned int m = reac->getListOfReactants()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            const std::string subId = reac->getReactant(j)->getSpecies();
            const std::string subComp = subId.substr(subId.size()-1, subId.size());
            const double subStoi = - reac->getReactant(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            initMatrix(reacCount, reacRev, subId, subComp, subStoi);
        }

//For loop for products
        m = reac->getListOfProducts()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            const std::string prodId = reac->getProduct(j)->getSpecies();
            const std::string prodComp = prodId.substr(prodId.size()-1, prodId.size());
            const double prodStoi = + reac->getProduct(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            initMatrix(reacCount, reacRev, prodId, prodComp, prodStoi);
        }
//Store name and index for reverse reactions
        if(reacRev) {
            reacName[fwdCount + (revCount - 1)] = reacId + "_reverse";
            ++revCount;
        }
    }
}

void Bacterium::initMatrix(const unsigned int &reacCount, const bool &reacRev,
                             const std::string &Id, const std::string &Comp, const double &Stoi) {
//##############################//
//  GLPK MATRIX INITIALISATION  //
//##############################//
    if(Comp == "e") {
//Check if the compartment is environment
        const int Num = environment.find(Id)->second;
        vecMeta.push_back(Num + 1), vecReac.push_back(reacCount + 1), vecStoi.push_back(+Stoi);
        if(reacRev) {
            //check if reaction is reversible, if so create the reverse reaction
            vecMeta.push_back(Num + 1), vecReac.push_back(fwdCount + revCount), vecStoi.push_back(-Stoi);
        }
    }
    else if(Comp == "c") {
//Check if the compartment is cytosol
        const int Num = cytosol.find(Id)->second;
        vecMeta.push_back(Num + 1), vecReac.push_back(reacCount + 1), vecStoi.push_back(+Stoi);
        if(reacRev) {
            //check if reaction is reversible, if so create the reverse reaction
            vecMeta.push_back(Num + 1), vecReac.push_back(fwdCount + revCount), vecStoi.push_back(-Stoi);
        }
    }
    else {
//Otherwise throw an error
        std::string errorMessage = std::string("unexpected compartment whilst reading SMBL file for species: ") + speciesID;
        throw std::logic_error(errorMessage);
    }
}

void Bacterium::doFBA() {
//##############################//
// GLPK READ AND SOLVE FUNCTION //
//##############################//

//Create new GLPK problem
            glp_prob *lp;
            lp = glp_create_prob();
            glp_set_obj_dir(lp, GLP_MAX);

//Add row and columns for stoichiometry matrix
            glp_add_cols(lp, reacName.size());
            glp_add_rows(lp, metaName.size());

//Assign metabolite names and bounds
            for(unsigned int i = 0; i < metaName.size(); ++i) {
                glp_set_row_name(lp, i+1, metaName.find(i)->second.c_str());
                glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
            }

//Assign reaction names and bounds

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

//Load and solve GLPK problem with Simplex method
            //vector.data() to access the internal arrays
            //vector.data()-1 to skip the first[0] element
            glp_load_matrix(lp, vecStoi.size(), vecMeta.data()-1, vecReac.data()-1, vecStoi.data()-1);
            glp_simplex(lp, NULL);

//Delete the current GLPK problem
            glp_delete_prob(lp);
            std::cout << "Succesfully finished Linear Programming step.\n\n";
}

void Bacterium::outputMatrix() {
//##############################//
//    TEST OUTPUT FUNCTIONS     //
//##############################//

        std::cout << "Number of external metabolites: " << environment.size() << "\n";
        std::cout << "Number of internal metabolites: " << cytosol.size() << "\n";
        std::cout << "Total number metabolites: " << metaName.size() << "\n";
        std::cout << "\nNumber of forward reactions: " << fwdCount << "\n";
        std::cout << "Number of reverse reactions: " << revCount - 1 << "\n";
        std::cout << "Total number of reactions: " << fwdCount + (revCount - 1) << "\n";
        std::cout << "\nNumber of non-zeroes in sparse matrix: " << vecStoi.size() << "\n\n";

        //Output sparse matrix with reaction and metabolite IDs to .csv file
        std::string outputFile = "./Sparse_Matrices/SM_" + speciesID + ".csv";
//Check if Sparse_Matrices/ folder exists, otherwise make it
        std::ofstream matrix(outputFile);
        if(!matrix.is_open()) {
            std::string errorMessage = std::string("could not open output file: ") + outputFile;
            throw std::runtime_error(errorMessage);
        }

        matrix << "Matrix Index\tMetabolite#\tMetaboliteID\tReaction#\tReactionID\tStoichiometry\n";
        for(unsigned int i = 1u; i <= vecStoi.size(); ++i){
            matrix << i << "\t" << vecMeta[i-1] << "\t" << metaName.find(vecMeta[i-1]-1)->second
                        << "\t" << vecReac[i-1] << "\t" << reacName.find(vecReac[i-1]-1)->second
                        << "\t" << vecStoi[i-1] << "\n";
        }
        matrix.close();
        std::cout << "Succesfully created sparse matrix in file: " << outputFile << "\n\n";
}

#endif // BACTERIUM_H
