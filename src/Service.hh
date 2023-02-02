#include <Service.hh>
#include <map>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <any>
#include<charconv>

template<int index,typename T,typename... args>
struct TemplateLoop{
    template<void givenFunction(T* val,args... a)>
    inline static const void runGivenFunction(T* val,args... arguments){
        givenFunction<index>(val+1,(arguments,...));
        TemplateLoop<index-1,T,args...>::runGivenFunction<givenFunction>(val,(arguments,...));
    }
};

template<typename T,typename... args>
struct TemplateLoop<0,T,args...>{
    template<void givenFunction(T* val,args... arguments)>
    inline static const void runGivenFunction(T* val,args... arguments){
        
    }
};

template<typename T1,template<typename T2> class stencil>
struct GenerateMRT{
    static const int m_Q=stencil<T1>::m_Q;
    int m_MRTMatrix[m_Q][m_Q];
    
    constexpr GenerateMRT():m_MRTMatrix(){
        for (int i=0;i<m_Q;i++){
            setMRTMatrixColumn<0>(i);
        }
    }

    template<int columnidx>
    inline constexpr void setMRTMatrixColumn(int rowidx){
        m_MRTMatrix[rowidx][columnidx]=stencil<T1>::template m_Moments<columnidx>[rowidx];
        setMRTMatrixColumn<columnidx+1>(rowidx);
    }

    template<>
    inline constexpr void setMRTMatrixColumn<m_Q>(int i){
        
    }
};

struct MRT{

};

template <typename T>
struct Counter
{
    Counter() {++counter;}
    virtual ~Counter() {--counter;}
    static int counter;
};
template <typename classtocount> int Counter<classtocount>::counter(0);

//===========================================================
  // Will read a given input file and assign variable values
  // based on formatting specified previously
//===========================================================
void readInput(const std::string fileName);



//===========================================================
  // Function will set variables based on an input vector of
  // strings consisting of variable names and values.
//===========================================================
template<typename T>
inline void SetVar(const std::string input,const std::string val);



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
std::vector<std::string>stringToVar(std::string init_string);


void collide(float* one,int a){}
int in;
float *test1;

void test(){
    TemplateLoop<19,float,int>::runGivenFunction<collide>(test1,in);
}
/*
template<typename T, typename stencil>
template<int idx>
inline constexpr void GenerateMRT<T,stencil>::setM(int i){
        m_M[i][idx]=stencil::template m_Moments<idx>[i];
        setM<idx+1>(i);
}
*/

//===========================================================
  // Will read a given input file and assign variable values
  // based on formatting specified previously
//===========================================================
void readInput(const std::string fileName)
{
	std::ifstream inputFile;
  inputFile.open(fileName);                                          // Opens input file

  std::string segment;                                               // Stores a segment of the input file (in this case a whole line)
  std::vector<std::string> lineSplit;                                // Stores an array of lines split based on some criteria (which will be used to determine variable names and variable values)
  std::vector<std::string> seglist;                                  // List of all the segments

  while(std::getline(inputFile, segment))                            // Read each line in inputFile as segment
    {
      lineSplit=stringToVar(segment);                                // Split the line using stringToVar function
      seglist.insert(std::end(seglist),std::begin(lineSplit),        // Insert the split line into the seglist vector
                     std::end(lineSplit));
    }
  
  inputFile.close();                                                 // Close the input file
  for(int i=0; i < input.size(); i++){
    if (i%2==0){
      SetVar<int>(seglist[i],seglist[i+1]);                                  // Function which sets variables based on the ordering in the seglist vector
      SetVar<float>(seglist[i],seglist[i+1]);
      SetVar<double>(seglist[i],seglist[i+1]);
      SetVar<string>(seglist[i],seglist[i+1]);
      SetVar<long>(seglist[i],seglist[i+1]);
    }
  }  
}



//===========================================================
  // Function will set variables based on an input vector of
  // strings consisting of variable names and values.
//===========================================================
template<typename T>
inline void SetVar(const std::string input,const std::string val){
  if(input.size()%2!=0){
    std::cout<<"WARNING: Input file may not be constructed "
               "correctly. Ensure Variables are defined as: "
               "name=value and everything else is commented "
               "with \% or \#.\n\nYou can separate variables"
               " using spaces, tabs or new lines. The code "
               "is pretty lenient with spaces and tabs (will"
               " be ignored if between an equals sign and a " 
               "name/value, will be ignored if repeated "
               "etc)."<<std::endl;
  }
  std::map<std::string, T> Keys=simulation_parameter<T>::VarKeys;
  
  if ( Keys.find(input) == Keys.end()) {
      std::cout<<"Cannot find a variable named "
              <<input<<". This could mean the "
              "names dont match but if you're seeing a "
              "lot of these errors then the structure "
              "of your file is probably not right."
              <<std::endl;
  } 

  else {
    T *ptr=&(Keys[input]);
    std::from_chars(val.data(), val.data() + val.size(), ptr)
    if (myPE==0)   std::cout<<input<<" = "
                              <<val<<std::endl;
  }
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
std::vector<std::string>stringToVar(std::string init_string){

  std::vector<std::string> final_StringVector;                       // Final vector containing variable names and variable values
  init_string=init_string.substr(0,init_string.find("%", 0));        // Remove comments (%) from string
  char current;                                                      // Store current character
  int eqPos=0;                                                       // Will be updated to store the position of '='
  bool lookForVariable=false;                                        // Approach switches if we are looking for a variable after an '=' character
  std::string tempString;                                            // Temporary string that will be reused
  std::string::size_type StringSize=init_string.size();              // Store string size

  for (std::string::size_type i = 0;                                 // Iterate through string
       i < StringSize; i++) {

    current=init_string[i];                                          // Current character
    if(lookForVariable==false){                                      // False by default as we havent reached an '='
      
      while ((current==' '||current=='\t')&&i<StringSize){           // Skip spaces and tabs
        i++;
        current=init_string[i];
      }
      
      if (current=='='){                                             // '=' found. Add variable name to vector (text before '=', stored in tempString)
        final_StringVector.push_back(tempString);
        tempString.clear();
        lookForVariable=true;                                        // Change this bool to true so we move on to the else statement below
      }
      
      else{
        tempString+=current;                                         // Add characters (not '=') to temporary string
      }
    }

    else{                                                            // We now look for a variable after the '='
      while ((current==' '||current=='\t')&&i<StringSize){           // Ignore spaces and tabs
        i++;
        current=init_string[i];
      }
      
      if (current!=' '&&current!='\t'){                              // First non space character
        if (current=='"'){                                           // If it is a '"' character, we must treat it as a string until we encounter another '"'
          i++;                                                       // Skip so we don't store the '"'
          current=init_string[i];
          
          while(current!='"'&&i < StringSize){                       // Add everything inside '"' characters to tempString
            tempString+=current;
            i++;
            current=init_string[i];
          }
          i++;                                                       // Skip so we don't store the '"'
          current=init_string[i];
        }
        else{
          while((current!=' '&&current!='\t')&&i<StringSize){        // If not a string, add all charactrs to tempstring until there is a space
            tempString+=current;
            i++;
            current=init_string[i];
          }
        }
        final_StringVector.push_back(tempString);                    // Add tempString to string vector
        lookForVariable=false;                                       // Start looking for variable name again
        tempString.clear();                                          // Clear tempString
        
        while ((current==' '||current=='\t')&&i<StringSize){         // Skip spaces after variable
          i++;
          current=init_string[i];
        }
        i--;
      }
    }
  }
  return final_StringVector;                                         // Return vector of alternating variable names and values
}