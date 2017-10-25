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

    //Creating test output files
    std::ofstream output("test.txt");
    if(!output) {
        std::cerr << "Could not open output file." << std::endl;
        exit(1);
    }
    std::ofstream output2("test2.txt");
    if(!output2) {
        std::cerr << "Could not open output2 file." << std::endl;
        exit(1);
    }
    std::ofstream output3("test3.txt");
    if(!output3) {
        std::cerr << "Could not open output3 file." << std::endl;
        exit(1);
    }

    //Read test input file with SBML
    SBMLDocument* document = readSBML(filename);
    Model *model = document->getModel();
    unsigned int errors = document->getNumErrors();
    document->printErrors(std::cerr);

    ListOfSpecies *spList = model->getListOfSpecies();
    ListOfReactions *reacList = model->getListOfReactions();

    //Write basic information to output
    output << "filename: "   << filename << "\n";
    output << "validation error(s): "   << errors   << "\n";
    output << "number of metabolites:\t" << spList->size() << "\n";
    output << "number of reactions:\t" << reacList->size() << "\n";

    //Create a GLPK problem with size #reactions and #metabolites
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_cols(lp, reacList->size());
    glp_add_rows(lp, spList->size());

    //Create internal and external maps to link metabolite name to metabolite number
    std::map<std::string, int> environment;
    int envCount = 0;
    std::map<std::string, int> cytosol;
    int cytCount = 0;

    //Write metabolite and compartment data to output file
    output << "\nMETABOLITES\n\n";
    for(unsigned int i = 0u; i < spList->size(); ++i) {

        std::string spId = spList->get(i)->getId();
        std::string spName = spId.substr(0, spId.size()-2);
        std::string spComp = spId.substr(spId.size()-1, spId.size());
        const char *spName_char = spName.c_str();
        glp_set_row_name(lp, i+1, spName_char);

        //Store metabolite name and number in map
        if(spComp == "e") {
            environment[spName] = envCount;
            ++envCount;
        }
        else if(spComp == "c") {
            cytosol[spName] = cytCount;
            ++cytCount;
        }
        else {
            std::cerr << "Unexpected compartment whilst reading SBML file: " << filename << "\n";
            exit(2);
        }
        output << i << " : " << spName << " " << spComp << "\n";
    }

    //Write metabolite maps with key and values to ouput file 2
    output2 << "internal metabolites:\n";
    for(const auto &c : cytosol){
        output2 << c.first << " : " << c.second << "\n";
    }
    output2 << "\nexternal metabolites:\n";
    for(const auto &e : environment){
        output2 << e.first << " : " << e.second << "\n";
    }

    //Create arrays for GLPK stoichiometry matrix
    int ia[1+10000], ja[1+10000]; //need a good way to determine the size of the array from metabolite and reaction data...?
    double ar[1+10000];
    int matrixCount = 1;

    //Write reactions to output file
    output << "\nREACTIONS\n\n";
    for(unsigned int i = 0u; i < reacList->size(); ++i) {
        Reaction *reac = reacList->get(i);

        std::string reacId = reac->getId();
        const char *reacId_char = reacId.c_str();
        glp_set_col_name(lp, i+1, reacId_char);

        output3 << i << " : " << reacId << ":\n";

        output << "reaction " << i << " : " << reacId
               << " (" << reac->getName() << ")\n"
               << "reversible? "<< (reac->getReversible() ? "Y\n" : "N\n")
               << "fast? " << (reac->getFast() ? "Y\n" : "N\n");

        //For loop for substrates
        unsigned int m = reac->getListOfReactants()->size();
        output << "Reactants: " << m << "\n";
        for(unsigned int j = 0u; j < m; ++j) {
            std::string subId = reac->getReactant(j)->getSpecies();
            std::string subName = subId.substr(0,subId.size()-2);
            std::string subComp = subId.substr(subId.size()-1, subId.size());
            int subStoi = reac->getReactant(j)->getStoichiometry();
            output << (j == 0 ? "" : " + ");
            output   << reac->getReactant(j)->getStoichiometry()
                     << " " << subName << " " << subComp;

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(subComp == "e") {
                auto subNum = environment.find(subName);
                ia[matrixCount] = subNum->second+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
                output3 << subComp << subNum->first << ", " << subNum->second << " ";
            }
            else if(subComp == "c") {
                auto subNum = cytosol.find(subName);
                ia[matrixCount] = subNum->second+1, ja[matrixCount] = i+1, ar[matrixCount] = -subStoi;
                matrixCount++;
                output3 << subComp << subNum->first << ", " << subNum->second << " ";
            }
            else {
                std::cerr << "Unexpected compartment whilst reading reactants from SBML file: " << filename << "\n";
                exit(3);
            }
        }
        output3 << "--> ";

        //For loop for products
        m = reac->getListOfProducts()->size();
        output << "\nProducts: " << m << "\n";
        for(unsigned int j = 0u; j < m; ++j) {
            std::string prodId = reac->getProduct(j)->getSpecies();
            std::string prodName = prodId.substr(0, prodId.size()-2);
            std::string prodComp = prodId.substr(prodId.size()-1, prodId.size());
            int prodStoi = reac->getProduct(j)->getStoichiometry();
            output << (j == 0 ? "" : " + ");
            output   << reac->getProduct(j)->getStoichiometry()
                     << " " << prodName << " " << prodComp;

            //Build stoichiometry matrix using maps.find() to get metabolite number
            if(prodComp == "e") {
                auto prodNum = environment.find(prodName);
                ia[matrixCount] = prodNum->second+1, ja[matrixCount] = i+1, ar[matrixCount] = prodStoi;
                matrixCount++;
                output3 << prodComp << prodNum->first << ", " << prodNum->second << " ";
            }
            else if(prodComp == "c") {
                auto prodNum = cytosol.find(prodName);
                ia[matrixCount] = prodNum->second+1, ja[matrixCount] = i+1, ar[matrixCount] = prodStoi;
                matrixCount++;
                output3 << prodComp << prodNum->first << ", " << prodNum->second << " ";
            }
            else {
                std::cerr << "Unexpected compartment whilst reading products from SBML file: " << filename << "\n";
                exit(4);
            }

            //Doesn't read FBC gene products
            //Double list of products possible (normal products + gene products)
            //Moreover FBC contains info on flux bounds, and a flux objective
        }
        output3 << "\n\n";
        output << "\n\n";
    }

    //Cout sparse matrix
    for(int i=1; i < matrixCount; i++){
        std::cout << ia[i] << ", " << ja[i] << ", " << ar[i] << "\n"; //last reaction not sparse...?
    }

    //Load and solve GLPK matrix
    glp_load_matrix(lp, matrixCount-1,ia,ja,ar);
    glp_simplex(lp, NULL);

    //Close output files
    output.close();
    output2.close();
    output3.close();

    //Delete SBML and GLPK elements
    delete document;
    glp_delete_prob(lp);
    return 0;
}
