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

    //Create internal and external maps to link metabolite name to metabolite number
    std::map<std::string, int> environment;
    std::map<std::string, int> cytosol;
    unsigned int metaCount = 0u;

    //Read metabolite data from SBML file
    for(unsigned int i = 0u; i < spList->size(); ++i) {
        std::string spId = spList->get(i)->getId();
        std::string spName = spId.substr(0, spId.size()-2);
        std::string spComp = spId.substr(spId.size()-1, spId.size());
//        const char *spId_char = spId.c_str();
//        glp_set_row_name(lp, i+1, spId_char);
//        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);

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
    std::vector<std::pair<int, std::string>> vecMeta, vecReac;
    std::vector<double> vecStoi;

    std::pair<int, std::string> test(1, "test");
    std::pair<int, std::string> test2(2, "test2");

    vecMeta.push_back(test);
    vecMeta.push_back(test2);

    std::cout << vecMeta[0].first << "\t" << vecMeta[0].second;
    std::cout << (vecMeta.data())->second;

    //These counters are initialised at 1, because glp can't handle 0's
    unsigned int matrixCount = 1u;
    unsigned int revCount = 1u;

    //Write substrates and products from SMBL files into GLPK arrays
    for(unsigned int i = 0u; i < reacList->size(); ++i) {
        const Reaction *reac = reacList->get(i);
        const std::string reacId = reac->getId();
//        const char *reacId_char = reacId.c_str();
//        glp_set_col_name(lp, i+1, reacId_char);
//        glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, 1000.0);
//        glp_set_obj_coef(lp, i+1, 0.0);

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
//                vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
//                    vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
                    ++matrixCount;
                }
            }
            else if(subComp == "c") {
                const int subNum = cytosol.find(subName)->second;
//                vecMeta.push_back(subNum+1), vecReac.push_back(i+1), vecStoi.push_back(-subStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
//                    vecMeta.push_back(subNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(+subStoi);
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
//                vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
                    //create output file with reaction and metabolite name and numbers here, to analyse the internal/external reactions
//                    vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                    ++matrixCount;
                }
            }
            else if(prodComp == "c") {
//                if(prodName == "M_biomass") glp_set_obj_coef(lp, i+1, 1.0);
                const int prodNum = cytosol.find(prodName)->second;
//                vecMeta.push_back(prodNum+1), vecReac.push_back(i+1), vecStoi.push_back(+prodStoi);
                ++matrixCount;
                if(reac->getReversible()) { //check if reaction is reversible, if so create the reverse reaction
//                    vecMeta.push_back(prodNum+1), vecReac.push_back(reacList->size() + revCount), vecStoi.push_back(-prodStoi);
                    ++matrixCount;
                }
            }
            else {
                std::cerr << "Unexpected compartment whilst reading products from SBML file: " << filename << "\n";
                exit(4);
            }
        }
        if(reac->getReversible()) { //check if reaction is reversible, if so give reverse reaction a name, obj_coef, and bounds
//            glp_set_obj_coef(lp, reacList->size() + revCount, 0.0);
//            const std::string revReacId = reacId + "_reverse";
//            const char *reacId_char = revReacId.c_str();
//            glp_set_col_name(lp, reacList->size() + revCount, reacId_char);
//            glp_set_col_bnds(lp, reacList->size() + revCount, GLP_DB, 0.0, 1000.0);
            ++revCount;
        }
    }

    //Output sparse matrix with reaction and metabolite IDs to .csv file
    std::ofstream output("test.csv");
    output << "Matrix Index\tMetabolite#\tMetaboliteID\tReaction#\tReactionID\tStoichiometry\n";
    for(unsigned int i=1; i < matrixCount; ++i){
//        output << i << "\t" << vecMeta[i-1] //<< "\t" << glp_get_row_name(lp, vecMeta[i-1])
//                    << "\t" << vecReac[i-1] //<< "\t" << glp_get_col_name(lp, vecReac[i-1])
//                    << "\t" << vecStoi[i-1] << "\n";
    }
    output.close();

//    //Create a GLPK problem with size #reactions and #metabolites
//    glp_prob *lp;
//    lp = glp_create_prob();
//    glp_set_obj_dir(lp, GLP_MAX);
//    glp_add_cols(lp, reacList->size()*2); //times 2 to account for potential reversible reactions
//    glp_add_rows(lp, spList->size());

 //    glp_del_cols(lp, (reacList->size()*2 - (reacList->size() + (revCount - 1)), const int delCol[]));
//    Please advise how to remove the columns...

    std::cout << "\nNumber of external metabolites: " << environment.size() << "\n";
    std::cout << "Number of internal metabolites: " << cytosol.size() << "\n";
    std::cout << "Total number metabolites: " << spList->size() << "\n";
    std::cout << "\nNumber of forward reactions: " << reacList->size() << "\n";
    std::cout << "Number of reverse reactions: " << revCount - 1 << "\n";
    std::cout << "Total number of reactions: " << reacList->size() + (revCount - 1) << "\n"; //this can be used to correctly set glp_add_col
    std::cout << "\nNumber of non-zeroes in sparse matrix: " << matrixCount - 1 << "\n\n";

    //Load and solve GLPK matrix
    //vector.data() to access the internal arrays
    //vector.data()-1 to skip the first[0] element
//    glp_load_matrix(lp, matrixCount-1,vecMeta.data()-1,vecReac.data()-1,vecStoi.data()-1);
//    glp_simplex(lp, NULL);

    //Delete SBML and GLPK elements
    delete document;
//    glp_delete_prob(lp);
    return 0;
}
