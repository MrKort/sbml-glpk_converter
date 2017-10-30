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
    glp_add_cols(lp, reacList->size());
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
    double ar[1+10000];
    int matrixCount = 1;
    int reverCount = 1;

    //Write substrates and products from SMBL files into GLPK arrays
    for(unsigned int i = 0u; i < reacList->size(); ++i) {
        Reaction *reac = reacList->get(i);
        std::string reacId = reac->getId();
        const char *reacId_char = reacId.c_str();
        glp_set_col_name(lp, i+1, reacId_char);
        glp_set_col_bnds(lp, i+1, GLP_DB, -1000.0, 1000.0); //Negative bounds were not allowed? But they seem to work.. Test with disected reversible reactions

        if(i+1 != reacList->size()) {
            glp_set_obj_coef(lp, i+1, 0.0);
        }
        else {
            glp_set_obj_coef(lp, i+1, 1.0);
        }

        //For loop for substrates
        unsigned int m = reac->getListOfReactants()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            std::string subId = reac->getReactant(j)->getSpecies();
            std::string subName = subId.substr(0,subId.size()-2);
            std::string subComp = subId.substr(subId.size()-1, subId.size());
            double subStoi = reac->getReactant(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(subComp == "e") {
                int subNum = environment.find(subName)->second;
                ia[matrixCount] = subNum+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
            }
            else if(subComp == "c") {
                int subNum = cytosol.find(subName)->second;
                ia[matrixCount] = subNum+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
            }
            else {
                std::cerr << "Unexpected compartment whilst reading reactants from SBML file: " << filename << "\n";
                exit(3);
            }
        }

        //For loop for products
        m = reac->getListOfProducts()->size();
        for(unsigned int j = 0u; j < m; ++j) {
            std::string prodId = reac->getProduct(j)->getSpecies();
            std::string prodName = prodId.substr(0, prodId.size()-2);
            std::string prodComp = prodId.substr(prodId.size()-1, prodId.size());
            double prodStoi = reac->getProduct(j)->getStoichiometry();

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(prodComp == "e") {
                int prodNum = environment.find(prodName)->second;
                ia[matrixCount] = prodNum+1, ja[matrixCount] = i+1, ar[matrixCount] = +prodStoi;
                matrixCount++;
            }
            else if(prodComp == "c") {
                int prodNum = cytosol.find(prodName)->second;
                ia[matrixCount] = prodNum+1, ja[matrixCount] = i+1, ar[matrixCount] = +prodStoi;
                matrixCount++;
            }
            else {
                std::cerr << "Unexpected compartment whilst reading products from SBML file: " << filename << "\n";
                exit(4);
            }
        }
    }

    //Cout sparse matrix
    for(int i=1; i < matrixCount; i++){
        std::cout << i << ", " << ia[i] << ", " << ja[i] << ", " << ar[i] << "\n";
    }
    std::cout << "\n";

    //Load and solve GLPK matrix
    glp_load_matrix(lp, matrixCount-1,ia,ja,ar); //rename ia[],ja[],ar[] arrays to logical names
    glp_simplex(lp, NULL);

    //Delete SBML and GLPK elements
    delete document;
    glp_delete_prob(lp);
    return 0;
}
