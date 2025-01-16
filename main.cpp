#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <cctype>
#include <sstream>
#include <format>
#include <regex>
#include <getopt.h>
#include "pugixml.hpp"

using string_list = std::vector<std::vector<std::string>>;

std::vector<std::string> string_to_vector(const std::string& AEs){
    std::vector<std::string> result;
    std::stringstream ss(AEs);
    std::string token;

    while (std::getline(ss, token, ',')) {
        result.push_back(token);
    }

    return result;
}

std::set<int> string_to_set(const std::string& AEs){
    std::set<int> result;
    std::stringstream ss(AEs);
    std::string token;

    while (std::getline(ss, token, ',')) {
        result.insert(std::stoi(token));
    }

    return result;
}

class patient{
public:
    patient(const std::string& id, const std::vector<std::string>& substance_list,
            const std::vector<std::string>& AE_list);

    patient(const std::string& id, const std::string& code_ATC,
            const std::string& AE_list);

    inline void set_substance_list(const std::vector<std::string>& substance_list){
        substance_list_ = substance_list;
    };

    inline void set_AE_list(const std::vector<std::string>& AE_list){
        AE_list_ = AE_list;
    };

    inline void set_ATC_code_list(const std::set<int>& ATC_code_list){
        ATC_code_list_ = ATC_code_list;
    };

    inline std::vector<std::string> get_substance_list() const{
        return substance_list_;
    }

    inline std::vector<std::string> get_AE_list() const{
        return AE_list_;
    }

    inline std::set<int> get_code_list() const{
        return ATC_code_list_;
    }

    inline void print_code(std::ostream& ofs) const {
        std::string result;
        for (auto it = ATC_code_list_.begin(); it != ATC_code_list_.end(); ++it) {
            result += std::to_string(*it);
            if (std::next(it) != ATC_code_list_.end()) {
               result += ':';
            }
        }
        ofs << result;
    }

    inline void print_AE(std::ostream& ofs) const{
        for(int i = 0 ; i < AE_list_.size()-1; ++i)
            ofs << AE_list_[i] << ",";
        ofs << AE_list_[AE_list_.size()-1];
            
    }

    inline void print_substance(std::ostream& ofs) const{
        for(auto& sub : substance_list_)
            ofs << sub << ' ';
    }

    void export_csv(std::ofstream& ofs) const;

private:
    std::string id_;
    std::vector<std::string> substance_list_;
    std::vector<std::string> AE_list_;
    std::set<int> ATC_code_list_;
};

patient::patient(const std::string& id, const std::vector<std::string>& substance_list,
            const std::vector<std::string>& AE_list) : id_{id}, substance_list_{substance_list},
            AE_list_{AE_list}
            {}

patient::patient(const std::string& id, const std::string& code_ATC,
            const std::string& AE_list) : id_{id}, ATC_code_list_{string_to_set(code_ATC)},
            AE_list_{string_to_vector(AE_list)}
            {}

void patient::export_csv(std::ofstream& ofs) const{
    print_code(ofs);
    ofs << ";";
    print_AE(ofs);
    ofs << ";";
    print_substance(ofs);
    ofs<<";";
}


bool load_xml_file(const std::string& path, pugi::xml_document& doc){
    
    pugi::xml_parse_result result = doc.load_file(path.c_str());

    if(!result){
        std::cout << "Error parsing the xml file.\n";
        return false;
    }
    return true;

}

//for each patient return a drug set in the form drug1;drug2; ... 
std::map<std::string,std::string> patient_drugs(const pugi::xpath_node_set& patientNodes,
                                       const std::vector<std::string>& id_list){

    std::map<std::string,std::string> returned_map;
    
    int i = 0;
    for(pugi::xpath_node patient : patientNodes){

        //recuperer les médicaments pris par le patient
        auto drugs = patient.node().children("drug");
        std::string drug_set = "";
        std::string d_name ="";
        
        std::string id = id_list[i++];
        for (const auto& drug : drugs){
            d_name = drug.child("medicinalproduct").child_value();
            std::transform(d_name.cbegin(), d_name.cend(), d_name.begin(),
                                            [](unsigned char c){ return std::tolower(c);});
            drug_set.append(d_name + ";");
        }
        returned_map.insert_or_assign(id,drug_set);
    }
    return returned_map;
}

//for each patient return an adverse event set in the form {AE1, AE2, ... , AEn}
std::map<std::string, std::vector<std::string>> patient_adverse_events(const pugi::xpath_node_set& patientNodes,
                                    const std::vector<std::string>& id_list){
    std::map<std::string, std::vector<std::string>> returned_set;
    
    std::vector<std::string> patient_AE;
    int i = 0;
    for(pugi::xpath_node patient : patientNodes){
        std::string AE_name = "";
        std::string id = id_list[i++];
        //catch all reaction for each patients
        auto AEs = patient.node().children("reaction");
        for (const auto& AE : AEs){
            //here we may want to lower case the PT term ? can do it post
            std::string tm = AE.child("reactionmeddrapt").child_value();
            std::transform(tm.begin(), tm.end(), tm.begin(),
                    [](unsigned char c){ return std::tolower(c); });
            
            patient_AE.push_back(tm);            
        }
        returned_set.insert_or_assign(id,patient_AE);
        patient_AE.clear();
    }
    return returned_set;
}


void write_tree_csv(std::ofstream& ost, const std::map<std::string,int>& tree){
    if(!ost.is_open()){
        std::cout<< "Error opening patient CSV file.\n";
        return; 
    }
    ost<< "ATCCode,Name,ATC_length\n";
}

std::map<std::string,std::string> get_standardized_substance(std::ifstream& ist){
    std::map<std::string,std::string> returned_map;
    if(!ist.is_open()){
        std::cout << "Error opening drug_standardized.csv\n";
        return returned_map;
    }

    std::string delimiter = ";";
    std::string quote = "\"";
    std::string drug = "";
    std::string substances = "";
    int pos;
    for(std::string line; std::getline(ist, line); ){
        pos = line.find(delimiter);
        drug = line.substr(0, pos);
        line.erase(0,pos+delimiter.length());
        //when there is multiple substances for a drug there are in the form : "sub1;sub2..." and we dont want the quotes to be here
        if((pos = line.find(quote)) != std::string::npos){
            line = line.substr(pos + quote.length(), line.find(quote,pos+quote.length())-1);
        }
        returned_map.insert({drug,line});
    }

   return returned_map;

}

std::map<std::string,std::string> get_atc_from_standardized(std::string_view path){
    std::map<std::string,std::string> returned_map;
    
    std::ifstream ist(path);
    if(!ist.is_open()){
        std::cerr << "Error opening the ATC binder file.\n";
        return returned_map;
    }

    std::string header;
    std::getline(ist, header);

    std::string delimiter = ";";
    int pos;
    std::string substance;
    std::vector<std::string> current_row;

    for(std::string line; std::getline(ist, line); ){
        while ((pos = line.find(delimiter)) != std::string::npos) {
            current_row.push_back(line.substr(0, pos));
            line.erase(0, pos + delimiter.length());
        }
        //we are interested in the susbtance name + ATC primary code so index 1 and 3 of the csv file
        returned_map.insert({current_row[1],current_row[3]});
        current_row.clear();
    }


    return returned_map;
}

std::map<std::string,std::vector<std::string>> extract_drugs_from_raw(const std::map<std::string,std::string>& raw_drugs){
    std::string delim = ";";
    std::map<std::string,std::vector<std::string>> patients_drugs_list;

    std::string current_drug = "";
    std::vector<std::string> drugs_list;
    for(auto [k,v] : raw_drugs){
        size_t pos = 0;
        //separe each drug
        while((pos = v.find(delim)) != std::string::npos){
            current_drug = v.substr(0,pos);
            v.erase(0, pos + delim.length());
            drugs_list.push_back(current_drug);
        }
        patients_drugs_list.emplace(k,drugs_list);
        drugs_list.clear();
    }
    return patients_drugs_list;
}

//apply a correction to a drug, in order to find a match in the drug-substances dictionnary
//we lemmatize the drug in parameter.
std::string apply_correction_drug(const std::string& drug){
    // we should be able to catch ~90% of uncorrect words
    std::string corrected_drug = drug;
    //if the delimiter is in the string, correct it
    if(drug.find_first_of("/([{^") != std::string::npos){
        corrected_drug = drug.substr(0,drug.find(" "));
    }
    return corrected_drug;
}

//some substances have an extra space at the end of the substances which 
//falsify every comparaison with atc tree, 1st idea -> erase this space
//2nd idea, some substances are represented with '-', int the ATC tree 
//there are no - , we then remove it
std::string apply_correction_substances(const std::string& substance){
    std::string corrected_substance = substance;
    if(substance.ends_with("\r") ){
        corrected_substance = substance.substr(0,substance.length()-1);
    }
    
    /*auto pos = corrected_substance.find("-");
    if(pos != std::string::npos){
        corrected_substance.erase(std::remove(corrected_substance.begin(),
                                              corrected_substance.end(),'-'),
                                  corrected_substance.end());        
    }*/
    //if(corrected_substance.ends_with(",") || corrected_substance.ends_with(","))
    return corrected_substance;
}

void get_patients_substances(
    const std::map<std::string,std::string>& standardized_dic,
    std::vector<patient>& patient_list)
    {
    string_list all_patients_substances;
    all_patients_substances.reserve(patient_list.size());
    
    std::string delimiter = ";";
    std::string current_drug = "";
    std::string substances = "";
    std::string current_sub;
    int err_catch = 0 ;
    std::vector<std::string> patient_substances;
    
    for(auto& patient : patient_list){
        //for each patients drugs
        for(auto drug : patient.get_substance_list()){
            try{
                //we apply basic drug correction in order to find a matching in the map
                substances = standardized_dic.at(apply_correction_drug(drug));
                substances = apply_correction_substances(substances); 
            }
            catch(const std::out_of_range& oor){
                substances = "NA";
            }
            int pos = substances.find(delimiter);
            while(pos != std::string::npos){
                current_sub = substances.substr(0,pos);
                substances.erase(0,pos + delimiter.length());
                patient_substances.push_back(current_sub);
                pos = substances.find(delimiter);
            }
            //add the last substances
            patient_substances.push_back(substances);
        }
        patient.set_substance_list(patient_substances);
        patient_substances.clear();
    }
}

std::vector<patient> delete_NA_substances(const std::vector<patient>& patients_list){
    std::vector<patient> returned_patients;
    returned_patients.reserve(patients_list.size());

    for(const auto& pat : patients_list){
        auto subs = pat.get_substance_list();
        if(std::count(subs.begin(), subs.end(), "NA") == 0){
            returned_patients.push_back(pat);
        }
    }
    returned_patients.shrink_to_fit();
    return returned_patients;
}

std::vector<patient> delete_NA_code(const std::vector<patient>& patients_list){
    std::vector<patient> returned_patients;
    returned_patients.reserve(patients_list.size());

    for(const auto& pat : patients_list){
        auto codes = pat.get_code_list();
        if(std::count(codes.begin(), codes.end(), INT_MIN) == 0){
            returned_patients.push_back(pat);
        }
    }
    returned_patients.shrink_to_fit();
    return returned_patients;
}

std::map<std::string, int> get_atc_code(std::ifstream& ist){
    std::map<std::string,int> atc_codes;
    if(!ist.is_open()){
        std::cout << "Error opening the ATC tree file : ATC_tree.csv \n";
        return atc_codes;
    }
    //get the header of the CSV file
    std::string header;
    std::getline(ist, header);

    std::string ATC_code;
    std::string sep = ",";
    int index = 0;
    int pos;
    //second columns of the csv file
    for(std::string line; std::getline(ist, line); ){
        pos = line.find(sep);
        ATC_code = line.substr(0,pos);
        atc_codes.insert({ATC_code, index++});
    }

    return atc_codes;
}
//convert substances to ATC code to allow data to be processed by the algorithm
std::vector<std::vector<int>> get_ATC_code_from_substances(const string_list& substances_lists,
                                                           const std::vector<std::string>& ATC_sub){
    std::vector<std::vector<int>> ATC_patients;
    ATC_patients.reserve(substances_lists.size());

    std::vector<int> patient;
    for(const auto& list : substances_lists){
        for(const auto& substances : list){
            //get the index of each substances i the ATC tree
            auto i = std::find(ATC_sub.begin(), ATC_sub.end(), substances);
            //if the substance has not been found, we add INT_MIN to mark the substances
            int index = i == ATC_sub.end() ? INT_MIN : std::distance(std::begin(ATC_sub), i);
            patient.push_back(index);
        }
        ATC_patients.push_back(patient);
        patient.clear();
    }
    return ATC_patients;
}

//convert substances to ATC code to allow data to be processed by the algorithm
void get_ATC_code_from_substances(std::vector<patient>& patients,
                                  const std::map<std::string,std::string>& map_ATC){
    std::vector<patient> ATC_patients;
    ATC_patients.reserve(patients.size());

    std::vector<std::string> patient;
    for(auto& pat : patients){
        for(const auto& substances : pat.get_substance_list()){
            try{
                patient.push_back(map_ATC.at(substances));
            }catch(const std::out_of_range& oor){
                patient.push_back("NA");
            }
        }
        pat.set_substance_list(patient);
        patient.clear();
    }
}

void get_index_from_ATC_code(std::vector<patient>& patients,
                                                      const std::map<std::string,uint16_t>& map_ATC_index){
    std::vector<patient> returned_pat;
    returned_pat.reserve(patients.size());

    std::set<int> patient_code;
    int cd;

    for(auto& patient : patients){
        for(const auto& code : patient.get_substance_list()){
            try{
                cd = map_ATC_index.at(code);
            }catch(const std::out_of_range& oor){
                cd = INT_MIN;
            }
            patient_code.insert(cd);
        }
        patient.set_ATC_code_list(patient_code);
        patient_code.clear();
    }
}

std::map<std::string,std::string> get_atc_tree(std::ifstream& ist){
    std::map<std::string,std::string> atc_line;
    if(!ist.is_open()){
        std::cout << "Error opening the ATC file ATC_tree.csv\n";
        return atc_line;
    }
    //get the header of the CSV file
    std::string header;
    std::getline(ist, header);

    std::string ATC_code, ATC_name;
    std::string sep = ",";
    int index = 0;
    int pos;
    //second columns of the csv file
    for(std::string line; std::getline(ist, line); ){
        pos = line.find(sep);
        ATC_code = line.substr(0,pos);
        ATC_name = line.substr(pos+sep.length(), line.size()-1);
        atc_line.insert({ATC_code, ATC_name});
    }

    return atc_line;
}
std::map<std::string, uint32_t> get_number_of_AE_per_AE(const string_list& patients_AE){
    std::map<std::string, uint32_t> AE_count;

    for(const auto& vec : patients_AE){
        for(const std::string& AE : vec){
            if(AE_count.contains(AE)){
                AE_count.at(AE)++;
            }else{
                AE_count.insert({AE,1});
            }
        }
    }

    return AE_count;
}

std::map<std::string, uint16_t> get_atc_tree_index(std::ifstream& ist){
    std::map<std::string, uint16_t> atc_line;
    uint16_t ATC_index= 0;
    if(!ist.is_open()){
        std::cout << "Error opening the ATC tree file ATC_tree.csv\n";
        return atc_line;
    }

    //get the header of the CSV file
    std::string header;
    std::getline(ist, header);

    std::string ATC_code;
    std::string sep = ",";
    int pos;
    //second columns of the csv file
    for(std::string line; std::getline(ist, line); ){
        pos = line.find(sep);
        ATC_code = line.substr(0,pos);
        atc_line.insert({ATC_code, ATC_index++});
    }

    return atc_line;
}

void export_patients(const std::vector<patient>& clean_patients_list, std::string_view out_path){
    std::ofstream ofs(out_path);
    if(!ofs.is_open()){
        std::cout << "Error opening: " << out_path <<  "\n";
        return;
    }
    ofs << "CODE ; AE ; SUBSTANCES \n";
    for(const auto& pat : clean_patients_list){
        pat.export_csv(ofs);
        ofs<<'\n';
    }

    ofs.close();
}

void export_code_with_AE(const std::vector<patient>& clean_patients_list,
                             const std::vector<bool>& AE, std::string_view out_path){
    std::ofstream ofs(out_path);
    if(!ofs.is_open()){
        std::cout << "Error opening: " << out_path <<  "\n";
        return;
    }
    ofs << "patientATC ; patientADR \n";
    for(int i = 0 ; i < clean_patients_list.size(); ++i){
        clean_patients_list[i].print_code(ofs);
        if(AE[i])
            ofs << ";1\n";
        else
            ofs << ";0\n";
    }

    ofs.close();
}


std::vector<patient> read_patients_csv(std::string_view in_path){
    std::vector<patient> returned_pat;
    std::ifstream ist(in_path);
    if(!ist.is_open()){
        std::cerr << "Error opening the patients csv file: "<< in_path << "\n";
        return returned_pat;
    }

    std::string header;
    std::getline(ist, header);

    std::string delimiter = ";";
    int pos;
    std::string substance;
    std::vector<std::string> current_row;
    int i = 0;
    for(std::string line; std::getline(ist, line); ){
        while ((pos = line.find(delimiter)) != std::string::npos) {
            current_row.push_back(line.substr(0, pos));
            line.erase(0, pos + delimiter.length());
        }
        

        returned_pat.emplace_back(std::to_string(i++), current_row[0], current_row[1]);
        current_row.clear();
    }
    return returned_pat;
}

std::vector<bool> get_AE_boolean(const string_list& patients_PT_code, std::string_view desired_PT){
    std::vector<bool> AE_true;
    AE_true.reserve(patients_PT_code.size());
    std::string desired_PT_lower = desired_PT.data();
    std::transform(desired_PT_lower.begin(), desired_PT_lower.end(),
                     desired_PT_lower.begin(),
                    [](unsigned char c){ return std::tolower(c); });
    
    for(const auto& patients : patients_PT_code){
        AE_true.push_back(std::find(patients.begin(), patients.end(), desired_PT_lower) != patients.end()
                            ? true
                            : false);
    }
    return AE_true;
}

std::vector<bool> get_AE_boolean_regex(const string_list& patients_PT_code, const std::regex& desired_PT_regex){
    std::vector<bool> AE_true;
    AE_true.reserve(patients_PT_code.size());
    
    for(const auto& patients : patients_PT_code){
        int i = 0;
        while(i < patients.size() && !std::regex_search(patients[i], desired_PT_regex))
            ++i;

        AE_true.push_back(i<patients.size());
    }

    return AE_true;
}

string_list AE_string_list_from_patient_vector(const std::vector<patient>& patients){
    std::vector<std::vector<std::string>> returned_list;
    returned_list.reserve(patients.size());
    for(const auto& pat : patients){
        returned_list.push_back(pat.get_AE_list());
    }

    returned_list.shrink_to_fit();
    return returned_list;
}

std::regex build_regex(const std::string& AE){
    std::string regex_pattern = "^.*\\b" + AE + "\\b.*$";
    return std::regex(regex_pattern, std::regex_constants::icase);
}

//need to process some rows in the form drug;"substances1;...;substances_n";...
std::vector<std::string> parse_csv_line(const std::string& line, char delimiter) {
    std::vector<std::string> columns;
    std::string cell;
    bool insideQuotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char currentChar = line[i];

        if (currentChar == '"') {
            // toggle the insideQuotes flag
            insideQuotes = !insideQuotes;
            cell += currentChar;
        } else if (currentChar == delimiter && !insideQuotes) {
            // if the delimiter is outside quotes, it's a column break
            columns.push_back(cell);
            cell.clear();
        } else {
            // Append the character to the current cell
            cell += currentChar;
        }
    }
    // Add the last cell
    columns.push_back(cell);

    return columns;
}

//remove every columns from original standardized_drugnames except substances and drug
void process_csv_standardized_drugnames(const std::string& inputFileName, const std::string& outputFileName) {
    std::ifstream inputFile(inputFileName);
    std::ofstream outputFile(outputFileName);

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files." << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::vector<std::string> columns = parse_csv_line(line, ';');

        // Ensure there are enough columns in the line
        if (columns.size() < 3) {
            std::cerr << "Skipping a malformed line: " << line << std::endl;
            continue;
        }

        outputFile << columns[1] << ";" << columns[2] << "\n";
    }

    inputFile.close();
    outputFile.close();
}

int main(int argc, char* argv[]){

    struct option long_options[] = {
        {"all", no_argument, nullptr, 'a'},
        {"specific", required_argument, nullptr, 's'},
        {"csvspecific", required_argument, nullptr, 'c'},
        {"input",required_argument, nullptr, 'i'},
        {"output",required_argument, nullptr, 'o'},
        {"mapping", required_argument, nullptr, 'm'},
        {"verbose", no_argument, nullptr, 'v'},
        {nullptr,0,nullptr,0}
    };

    int opt;
    bool all, mapping_processed = false, verbose = false;
    std::string specific_AE;
    std::string csv_specific_AE;
    std::string input_file;
    std::string output_file;// csv_outputfile
    std::string mapping_path;
    std::vector<patient> clean_patients_list;
    while((opt = getopt_long(argc, argv, "aps:c:i:o:m:", long_options, nullptr)) != -1){
        switch (opt)
        {
        case 'a':
            all = true;
            break;
        case 'p':
            mapping_processed = true;
            break;
        case 's':
            specific_AE = optarg;
            break;
        case 'c':
            csv_specific_AE = optarg;
            break;
        case 'i':
            input_file = optarg;
            break;
        case 'o':
            output_file = optarg;
            break;
        case 'm':
            mapping_path = optarg;
            break;
        case 'v':
            verbose = true;
            break;
        case '?':
            std::cerr << "Unknown option or missing argument.\n";
            return 1;
        }
    }

    if (input_file.empty()) {
        std::cerr << "Error: The --input option is mandatory. Please add it.\n";
        return 1;
    }
    if (output_file.empty()) {
        std::cerr << "Error: The --output option is mandatory. Please add it.\n";
        return 1;
    }
    if((all || !specific_AE.empty()) && (mapping_path.empty() && !mapping_processed ) ){
        std::cerr << "Error: The --mapping option is mandatory when going from xml to csv. Please add it.\n";
        return 1;
    }

    int option_count = (all ? 1 : 0) + (!specific_AE.empty() ? 1 : 0) + (!csv_specific_AE.empty() ? 1 : 0);
    if (option_count > 1) {
        std::cerr << "Error: Only one of --all, --specific, or --csvspecific can be specified at a time.\n";
        return 1;
    }

      
    std::cout << "Input file: " << input_file << "\n";
    std::cout << "Output file: " << output_file << "\n";

    
    pugi::xml_document doc;

   if(all || !specific_AE.empty()){  
    auto status = load_xml_file(input_file,doc);
    if(!status)
        return -1;
    
    auto patients = doc.select_nodes("/ichicsr/safetyreport/patient");
    auto patients_id = doc.select_nodes("/ichicsr/safetyreport");
    
    std::vector<std::string> id_list;
    id_list.reserve(patients_id.size());
    for(auto pat : patients_id){
        std::string bop = pat.node().child("safetyreportid").child_value();
        id_list.push_back(bop);
    }
    

    std::map<std::string,std::string> raw_patient_drug_list;
    //Get the drugs for each patients in the format "drug1;drug2;...drug_n";
    raw_patient_drug_list = patient_drugs(patients, id_list);
 

    //Convert the raw patients list in a nice vector structure
    std::map<std::string,std::vector<std::string>> patients_drugs_list;

    //Get the drugs for each patient in the form ["drug1","drug2",...,"drug_n"]
    patients_drugs_list = extract_drugs_from_raw(raw_patient_drug_list);
    // Get the Adverse events for each patients :
    auto patients_AEs_list = patient_adverse_events(patients, id_list);



    //auto patients_AE_counts = get_number_of_AE_per_AE(patients_AEs_list);
    
    std::vector<patient> patients_list;
    patients_list.reserve(patients_drugs_list.size());
    for(const auto& [k,v] : patients_drugs_list){
        patients_list.emplace_back(k,v,patients_AEs_list.at(k));
    }
    
    
    //convert the standardized drugname csv in a cpp structure (map justifiée car à priori un drugname par substance)
    std::string mapping_path_2_columns = "./drugnames_standardized_2_columns.csv";
    if(!mapping_processed)
        process_csv_standardized_drugnames(mapping_path, mapping_path_2_columns);
    std::ifstream ist(mapping_path_2_columns);
    std::map<std::string,std::string> map_standardized = get_standardized_substance(ist);
    ist.close();


    //here we need to be cautious, there are multiple ATC code for a single substance
    // furthermore some primary code are Z (?)
    std::string path_atc_mapping = "./ATC_binder_2024.csv";
    std::map<std::string, std::string> map_standardized_to_ATC = get_atc_from_standardized(path_atc_mapping);

    
    std::string path_tree = "./ATC_tree.csv";    
    std::ifstream ist_tree(path_tree); 
    //auto ATC_code = get_atc_code(ist_tree);
    //auto ATC_code_name = get_atc_tree(ist_tree);
    auto ATC_code_index = get_atc_tree_index(ist_tree);
    ist_tree.close();
    
    //from the patient drug list, get the standardized substances names in order to use the ATCtree
    
    string_list patients_substances_list;
    //each patient have now a substances list corresponding to their medication,
    //we found the substances composing each medication in "drugnames_standardized_med_only.csv"
    get_patients_substances(map_standardized, patients_list);
    if(verbose)
        std::cout << "Patient number before cutting NA substance : " << patients_list.size() << '\n';
    //we just delete row when the substances hasn't been found -> not the best ? 
    clean_patients_list = delete_NA_substances(patients_list);

    if(verbose)
        std::cout << "Patient number after cutting NA substance : " << clean_patients_list.size() << '\n';
    
    
    get_ATC_code_from_substances(clean_patients_list, map_standardized_to_ATC);
    //remove patients containing "NA" in their ATC list 
    clean_patients_list = delete_NA_substances(clean_patients_list);
    if(verbose)
        std::cout << "Patient number after cutting NA ATC_code : " << clean_patients_list.size() << '\n';

    get_index_from_ATC_code(clean_patients_list, ATC_code_index);
    clean_patients_list = delete_NA_code(clean_patients_list);
    if(verbose)
        std::cout << "Patient number after removing INT_MIN from ATC_code : " << clean_patients_list.size() << '\n';

    //we export this if user wants the file with every AE
    if(all && !output_file.empty()){
        export_patients(clean_patients_list, output_file);
        std::cout << "Succesfully exported data to : "<< output_file <<"\n";
    }
    

   }


    if(!all){
        std::vector<patient> imported_patients;
        if(!csv_specific_AE.empty())
            imported_patients = read_patients_csv(input_file);
        else
            imported_patients = clean_patients_list;
        
        std::string desired_AE = csv_specific_AE.empty() ? specific_AE : csv_specific_AE;
        std::regex AE_reg = build_regex(desired_AE);

        std::vector<bool> ADR = get_AE_boolean_regex(AE_string_list_from_patient_vector(imported_patients), AE_reg);
        export_code_with_AE(imported_patients, ADR, output_file);
        std::cout << "Succesfully exported data to : "<< output_file <<"\n";
    }

    return 0; 

}