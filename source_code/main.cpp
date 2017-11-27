#include <glpk.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/fbcfwd.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
#include <fstream>
#include <map>
#include <regex>

LIBSBML_CPP_NAMESPACE_USE

int main (int argc, char* argv[])
{
    //Declare test input file
    const char* filename = "/home/kevin/Git/sbml-glpk_converter/test.xml";

//##############################//
//  SBML FILE READING FUNCTION  //
//##############################//

    //Read test input file with SBML
    SBMLDocument* document = readSBML(filename);
    Model *model = document->getModel();
    unsigned int errors = document->getNumErrors();
    std::cout << errors << " errors while reading file: " << filename << "\n";
    document->printErrors(std::cerr);

    ListOfSpecies *spList = model->getListOfSpecies();
    unsigned int numMeta = spList->size();
    ListOfReactions *reacList = model->getListOfReactions();
    unsigned int numReac = reacList->size();

    //Create internal and external maps to link metabolite name to metabolite number
    std::map<std::string, int> environment;
    std::map<std::string, int> cytosol;
    std::map<int, std::string> metaName;
    unsigned int metaCount = 0u;

    //Read metabolite data from SBML file
    for(unsigned int i = 0u; i < spList->size(); ++i) {
        std::string spId = spList->get(i)->getId();
        std::string spName = spId.substr(0, spId.size()-2);
        std::string spComp = spId.substr(spId.size()-1, spId.size());

        metaName[metaCount] = spId; //store metaNumber and metaId in a map independent of compartmentalisation

        //Store metabolite name and number in map
        if(spComp == "e") {
            environment[spName] = metaCount;
            ++metaCount;
        }
        else if(spComp == "c") {
            cytosol[spName] = metaCount;
            ++metaCount;
        }
        else {
            std::cerr << "Unexpected compartment whilst reading SBML file: " << filename << "\n";
            exit(1);
        }
    }
    if(!cytosol.find("M_biomass")->second) {
        std::cerr << "Error: could not find a biomass producing reaction.\n";
        exit(2);
    }

    //Create vectors for GLPK
    std::vector<int> vecMeta, vecReac; //vector with GLPK index value
    std::vector<double> vecStoi;
    std::map<int, std::string> reacName;

    //These counters are initialised at 1, because glp can't handle 0's
    unsigned int matrixCount = 1u;
    unsigned int revCount = 1u;
    unsigned int ObjFunc;

    //Write substrates and products from SMBL files into GLPK arrays
    for(unsigned int i = 0u; i < reacList->size(); ++i) {
        const Reaction *reac = reacList->get(i);
        const std::string reacId = reac->getId();
        reacName[i] = reacId;

        //Determine wheter the reaction name is biomass0, this is usually the objective function
        //In the future update the objective function determenation by using SBML::FBC:objective_function
        if(reac->getName() == "biomass0") {
            std::cout << "\nObjective function found!\nID: " << i << "\tName: " << reac->getId() << "\n\n";
            ObjFunc = i;
        }

        //For loop for substrates
        unsigned int m = reac->getListOfReactants()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            const std::string subId = reac->getReactant(j)->getSpecies();
            const std::string subName = subId.substr(0,subId.size()-2);
            const std::string subComp = subId.substr(subId.size()-1, subId.size());
            const double subStoi = reac->getReactant(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(subComp == "e") {
                const int subNum = environment.find(subName)->second;
                vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
                    ++matrixCount;
                }
            }
            else if(subComp == "c") {
                const int subNum = cytosol.find(subName)->second;
                vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
                    ++matrixCount;
                }
            }
            else {
                std::cerr << "Unexpected compartment whilst reading substrates from SBML file: " << filename << "\n";
                exit(3);
            }
        }

        //For loop for products
        m = reac->getListOfProducts()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            const std::string prodId = reac->getProduct(j)->getSpecies();
            const std::string prodName = prodId.substr(0, prodId.size()-2);
            const std::string prodComp = prodId.substr(prodId.size()-1, prodId.size());
            const double prodStoi = reac->getProduct(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(prodComp == "e") {
                const int prodNum = environment.find(prodName)->second;
                vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                    ++matrixCount;
                }
            }
            else if(prodComp == "c") {
                const int prodNum = cytosol.find(prodName)->second;
                vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                    ++matrixCount;
                }
            }
            else {
                std::cerr << "Unexpected compartment whilst reading products from SBML file: " << filename << "\n";
                exit(4);
            }
        }
        if(reac->getReversible()) { //check if reaction is reversible, if so give reverse reaction a name, obj_coef, and bounds
            reacName[numReac + revCount-1] = reacId + "_reverse";
            ++revCount;
        }
    }

    //Delete SBML and GLPK elements
//    delete cytosol;
//    delete environment;
    delete document;

//##############################//
// GLPK READ AND SOLVE FUNCTION //
//##############################//

//        CREATE A NEW GLPK PROBLEM
        glp_prob *lp;
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MAX);

//        ADD ROW AND COLUMNS FOR STOICHIOMETRY MATRIX
        glp_add_cols(lp, numReac + (revCount - 1));
        glp_add_rows(lp, numMeta);

//        ASSIGN METABOLITE NAME AND BOUNDS
        for(unsigned int i = 0; i < numMeta; ++i){
            glp_set_row_name(lp, i+1, metaName.find(i)->second.c_str());
            glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
        }

//        ASSIGN REACTION NAME AND BOUNDS
        for(unsigned int i = 0; i < numReac + (revCount - 1); ++i) {
            glp_set_col_name(lp, i+1, reacName.find(i)->second.c_str());
            glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, 1000.0);
            glp_set_obj_coef(lp, i+1, 0.0);

//        CHECK FOR FUCTION OBJECTIVE REACTION
            if(i == ObjFunc) {
                glp_set_obj_coef(lp, i+1, 1.0);
            }
        }

//        LOAD AND SOLVE GLPK PROBLEM
//        vector.data() to access the internal arrays
//        vector.data()-1 to skip the first[0] element
        glp_load_matrix(lp, matrixCount-1,vecMeta.data()-1,vecReac.data()-1,vecStoi.data()-1);
        glp_simplex(lp, NULL);

//        DELETE THE GLPK PROBLEM
        glp_delete_prob(lp);

//##############################//
//    TEST OUTPUT FUNCTIONS     //
//##############################//

    //Output sparse matrix with reaction and metabolite IDs to .csv file
    std::ofstream matrix("sparse_matrix.csv");
    matrix << "Matrix Index\tMetabolite#\tMetaboliteID\tReaction#\tReactionID\tStoichiometry\n";
    for(unsigned int i=1; i < matrixCount; ++i){
        matrix << i << "\t" << vecMeta[i-1] << "\t" << metaName.find(vecMeta[i-1]-1)->second
                    << "\t" << vecReac[i-1] << "\t" << reacName.find(vecReac[i-1]-1)->second
                    << "\t" << vecStoi[i-1] << "\n";
    }
    matrix.close();

    std::cout << "\nNumber of external metabolites: " << environment.size() << "\n";
    std::cout << "Number of internal metabolites: " << cytosol.size() << "\n";
    std::cout << "Total number metabolites: " << numMeta << "\n";
    std::cout << "\nNumber of forward reactions: " << numReac << "\n";
    std::cout << "Number of reverse reactions: " << revCount - 1 << "\n";
    std::cout << "Total number of reactions: " << numReac + (revCount - 1) << "\n";
    std::cout << "\nNumber of non-zeroes in sparse matrix: " << matrixCount - 1 << "\n\n";

    return 0;
}
