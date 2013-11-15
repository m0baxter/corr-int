#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "StringManipulators.hpp"


std::vector<std::string> split( const std::string &s, const std::string &sep ) {

   std::vector<std::string> tokens;
   
   size_t i = 0;
   size_t j = 0;
   
   //deal with first token:
   i = s.find(sep,i);
   tokens.push_back(s.substr(j,i-j));
   j = i;
   ++i;
   
   while ( j < std::string::npos ) {
   
      i = s.find(sep,i);
      tokens.push_back(s.substr(j+1,i-j-1));
      j = i;
      ++i;      
   }  

   //deal with repition of the token:
   tokens.erase( remove(tokens.begin(), tokens.end(), sep ), tokens.end());
   //remove empty strings created by multiple tokens:
   tokens.erase( remove(tokens.begin(), tokens.end(), "" ), tokens.end());
   
   return tokens;
}
