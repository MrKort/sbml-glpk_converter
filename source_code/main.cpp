#include <glpk.h>
#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/fbcfwd.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>
#include <fstream>
#include <map>
LIBSBML_CPP_NAMESPACE_USE

int main (int argc, char* argv[])
{
    //Declare test input file
    const char* filename   = "/home/kevin/Git/sbml-glpk_converter/test.xml";

    //Read test input file with SBML
    SBMLDocument* document = readSBML(filename);
    Model *model = document->getModel();
    unsigned int errors = document->getNumErrors();
    std::cout << errors << " errors while reading file: " << filename << "\n";
    document->printErrors(std::cerr);

    ListOfSpecies *spList = model->getListOfSpecies();
    ListOfReactions *reacList = model->getListOfReactions();

    //Create a GLPK problem with size #reactions and #metabolites
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_cols(lp, reacList->size()*2); //times 2 to account for potential reversible reactions
    glp_add_rows(lp, spList->size());

    //Create internal and external maps to link metabolite name to metabolite number
    std::map<std::string, int> environment;
    std::map<std::string, int> cytosol;
    int metaCount = 0;

    //Read metabolite data from SBML file
    for(unsigned int i = 0u; i < spList->size(); ++i) {
        std::string spId = spList->get(i)->getId();
        std::string spName = spId.substr(0, spId.size()-2);
        std::string spComp = spId.substr(spId.size()-1, spId.size());
        const char *spName_char = spName.c_str();
        glp_set_row_name(lp, i+1, spName_char);
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);

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
            exit(2);
        }
    }

    //Create arrays for GLPK
    int ia[1+10000], ja[1+10000]; //need a good way to determine the size of the array from metabolite and reaction data...?
    double ar[1+10000]; //after first run the exact sizes for these arrays are known, and can be used for further calculations
    int matrixCount = 1;
    int revCount = 1;

    //Write substrates and products from SMBL files into GLPK arrays
    for(unsigned int i = 0u; i < reacList->size(); ++i) {
        const Reaction *reac = reacList->get(i);
        const std::string reacId = reac->getId();
        const char *reacId_char = reacId.c_str();
        glp_set_col_name(lp, i+1, reacId_char);
        glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, 1000.0);

        if(i+1 != reacList->size()) {
            glp_set_obj_coef(lp, i+1, 0.0);
        }
        else {
            glp_set_obj_coef(lp, i+1, 1.0); //only the last reaction (biomass production) is set as an objective function
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
                ia[matrixCount] = subNum+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    ia[matrixCount] = subNum+1, ja[matrixCount] = reacList->size() + revCount, ar[matrixCount] = +subStoi;
                    matrixCount++;
                }
            }
            else if(subComp == "c") {
                const int subNum = cytosol.find(subName)->second;
                ia[matrixCount] = subNum+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    ia[matrixCount] = subNum+1, ja[matrixCount] = reacList->size() + revCount, ar[matrixCount] = +subStoi;
                    matrixCount++;
                }
            }
            else {
                std::cerr << "Unexpected compartment whilst reading reactants from SBML file: " << filename << "\n";
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
                ia[matrixCount] = prodNum+1, ja[matrixCount] = i+1, ar[matrixCount] = +prodStoi;
                matrixCount++;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    ia[matrixCount] = prodNum+1, ja[matrixCount] = reacList->size() + revCount, ar[matrixCount] = -prodStoi;
                    matrixCount++;
                }
            }
            else if(prodComp == "c") {
                const int prodNum = cytosol.find(prodName)->second;
                ia[matrixCount] = prodNum+1, ja[matrixCount] = i+1, ar[matrixCount] = +prodStoi;
                matrixCount++;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    ia[matrixCount] = prodNum+1, ja[matrixCount] = reacList->size() + revCount, ar[matrixCount] = -prodStoi;
                    matrixCount++;
                }
            }
            else {
                std::cerr << "Unexpected compartment whilst reading products from SBML file: " << filename << "\n";
                exit(4);
            }
        }
        if(reac->getReversible()) { //check if reaction is reversible, if so give reverse reaction a name, obj_coef, and bounds
            glp_set_obj_coef(lp, reacList->size() + revCount, 0.0);
            const std::string revReacId = reacId + "_reverse";
            const char *reacId_char = revReacId.c_str();
            glp_set_col_name(lp, reacList->size() + revCount, reacId_char);
            glp_set_col_bnds(lp, reacList->size() + revCount, GLP_DB, 0.0, 1000.0);
            revCount++;
        }
    }

    //Output sparse matrix to csv file
    std::ofstream output("test.csv");
    output << "Matrix Index\tMetabolite#\tReaction#\tStoichiometry\n";
    for(int i=1; i < matrixCount; i++){
        output << i << "\t" << ia[i] << "\t" << ja[i] << "\t" << ar[i] << "\n"; // this can be used to correctly set ia, ja, and ar sizes
    }
    output.close();

    std::cout << "Maximum number of reactions: " << reacList->size() + revCount - 1 << "\n\n"; //this can be used to correctly set glp_add_col

    //Load and solve GLPK matrix
    glp_load_matrix(lp, matrixCount-1,ia,ja,ar);
    glp_simplex(lp, NULL);

    //Delete SBML and GLPK elements
    delete document;
    glp_delete_prob(lp);
    return 0;
}
