#include <any>
#include <map>
#include <string>

#include "Mpi.hh"

class InputParameters {
   public:
    template <typename T>
    void addParameter(T& parameter, std::string name) {
        mVarTypes[name] = typeid(T).name();
        mVarKeys[name] = &parameter;
    }

    template <typename T>
    void addParameter(std::string name) {
        mVarTypes[name] = typeid(T).name();
        mVarKeys[name] = new T;  // This should be freed at the end
    }

    template <typename T>
    T& getParameter(std::string name) {
        return *std::any_cast<T*>(mVarKeys.at(name));
    }

    std::string& getParameterType(std::string name) { return mVarTypes.at(name); }

    template <typename T>
    const T& getParameter(std::string name) const {
        return *std::any_cast<const T*>(mVarKeys.at(name));
    }

    const std::string& getParameterType(std::string name) const { return mVarTypes.at(name); }

    void readInput(std::string fileName) {
        if (mpi.rank == 0) std::cout << "Reading parameters for simulation from " << fileName << "." << std::endl;

        std::ifstream inputFile;
        inputFile.open(fileName);  // Opens input file

        std::string segment;                 // Stores a segment of the input file (in this case a whole line)
        std::vector<std::string> lineSplit;  // Stores an array of lines split based on some criteria (which will be
                                             // used to determine variable names and variable values)
        std::vector<std::string> seglist;    // List of all the segments

        while (std::getline(inputFile, segment))  // Read each line in inputFile as segment
        {
            lineSplit = stringToVar(segment);                         // Split the line using stringToVar function
            seglist.insert(std::end(seglist), std::begin(lineSplit),  // Insert the split line into the seglist vector
                           std::end(lineSplit));
        }

        inputFile.close();  // Close the input file

        SetVar(mVarKeys, mVarTypes,
               seglist);  // Function which sets variables based on the ordering in the seglist vector
    }

    //===========================================================
    // Function will set variables based on an input vector of
    // strings consisting of variable names and values.
    //===========================================================
    void SetVar(std::map<std::string, std::any>& Keys, std::map<std::string, std::string>& Types,
                const std::vector<std::string>& input) {
        if (input.size() % 2 != 0) {
            std::cout << "WARNING: Input file may not be constructed correctly. Ensure Variables are defined as: "
                         "name=value and everything else is commented with % or #.\n\nYou can separate variables using "
                         "spaces, tabs or new lines. The code is pretty lenient with spaces and tabs (will be ignored "
                         "if between an equals sign and a name/value, will be ignored if repeated etc)."
                      << std::endl;
        }

        std::string tempinput;  // String that will temporarily store some portion of the vector

        for (int i = 0; i < (int)input.size(); i++) {
            if (i % 2 == 0) {  // Take only even (including zero) indices as these correspond to variable names.
                tempinput = input[i];

                if (Keys.find(tempinput) == Keys.end() &&
                    mpi.rank == 0) {  // If the variable name is NOT found in the Keys map (variable names are
                                      // automatically stored in this map if defined as a "simulation_parameter")
                    std::cout << "Cannot find a variable named " << tempinput
                              << ". This could mean the "
                                 "names dont match but if you're seeing a "
                                 "lot of these errors then the structure "
                                 "of your file is probably not right."
                              << std::endl;
                }
                // Else if the variable name is NOT found in the Types
                // map
                else if (Types.find(tempinput) == Types.end() &&
                         mpi.rank == 0) {  // If the variable type is NOT found in the Keys map (variable type are
                                           // automatically stored in this map if defined as a "simulation_parameter")
                    std::cout << "WARNING: Cannot find type of variable " << tempinput
                              << ". Check that names match"
                                 " in \"VarKeys()\" and \"VarTypes\" in "
                                 "initialise_simulation.cc."
                              << std::endl;

                } else {  // Variable name and type found
                          // Depending on the desired type, convert string to that type

                    const std::string TYPE = Types[tempinput];  // Extract the type info from the Types array

                    if (TYPE == "i") {  // Integer
                        *std::any_cast<int*>(Keys[tempinput]) = std::stoi(
                            input[i + 1]);  // Set that variable with the variable in the next entry in the vector
                    }

                    if (TYPE == "f") {  // float

                        *std::any_cast<float*>(Keys[tempinput]) = std::stof(input[i + 1]);
                    }

                    if (TYPE == "d") {  // double

                        *std::any_cast<double*>(Keys[tempinput]) = std::stod(input[i + 1]);
                    }

                    if (TYPE == "l") {  // Long
                        *std::any_cast<long*>(Keys[tempinput]) = std::stol(input[i + 1]);
                    }

                    if (TYPE ==
                        "NSt7__cxx1112basic_stringIcSt11char_traitsIc"  // String
                        "ESaIcEEE") {
                        *std::any_cast<std::string*>(Keys[tempinput]) = input[i + 1];
                    }

                    if (mpi.rank == 0) std::cout << tempinput << " = " << input[i + 1] << std::endl;
                }
            }
        }
        if (mpi.rank == 0) std::cout << std::endl;
    }

    //===========================================================
    // Function to read a string, will split into variable_name
    // variable, variable_name, variable etc.......
    // Will ignore anything after % character
    // Will look for variable after = character
    // Will repeat after each variable
    // This is somewhat convoluted to allow for loose
    // formatting andmultiple variable definitions on the same
    // line
    //===========================================================
    std::vector<std::string> stringToVar(std::string init_string) {
        std::vector<std::string> final_StringVector;  // Final vector containing variable names and variable values
        init_string = init_string.substr(0, init_string.find("%", 0));  // Remove comments (%) from string
        char current;                                                   // Store current character
        // int eqPos=0;                                                       // Will be updated to store the position
        // of '='
        bool lookForVariable = false;  // Approach switches if we are looking for a variable after an '=' character
        std::string tempString;        // Temporary string that will be reused
        std::string::size_type StringSize = init_string.size();  // Store string size

        for (std::string::size_type i = 0; i < StringSize; i++) {
            current = init_string[i];        // Current character
            if (lookForVariable == false) {  // False by default as we havent reached an '='

                while ((current == ' ' || current == '\t') && i < StringSize) {  // Skip spaces and tabs
                    i++;
                    current = init_string[i];
                }

                if (current == '=') {  // '=' found. Add variable name to vector (text before '=', stored in tempString)
                    final_StringVector.push_back(tempString);
                    tempString.clear();
                    lookForVariable = true;  // Change this bool to true so we move on to the else statement below
                }

                else {
                    tempString += current;  // Add characters (not '=') to temporary string
                }
            }

            else {  // We now look for a variable after the '='
                while ((current == ' ' || current == '\t') && i < StringSize) {  // Ignore spaces and tabs
                    i++;
                    current = init_string[i];
                }

                if (current != ' ' && current != '\t') {  // First non space character
                    if (current ==
                        '"') {  // If it is a '"' character, we must treat it as a string until we encounter another '"'
                        i++;    // Skip so we don't store the '"'
                        current = init_string[i];

                        while (current != '"' &&
                               i < StringSize) {  // Add everything inside '"' characters to tempString
                            tempString += current;
                            i++;
                            current = init_string[i];
                        }
                        i++;  // Skip so we don't store the '"'
                        current = init_string[i];
                    } else {
                        while ((current != ' ' && current != '\t') &&
                               i < StringSize) {  // If not a string, add all charactrs to tempstring until there is a
                                                  // space
                            tempString += current;
                            i++;
                            current = init_string[i];
                        }
                    }
                    final_StringVector.push_back(tempString);  // Add tempString to string vector
                    lookForVariable = false;                   // Start looking for variable name again
                    tempString.clear();                        // Clear tempString

                    while ((current == ' ' || current == '\t') && i < StringSize) {  // Skip spaces after variable
                        i++;
                        current = init_string[i];
                    }
                    i--;
                }
            }
        }
        return final_StringVector;  // Return vector of alternating variable names and values
    }

   private:
    std::map<std::string, std::any> mVarKeys;
    std::map<std::string, std::string> mVarTypes;
};
