#include <string>
#include <vector>

#ifndef STRINGMANIP_HPP
#define STRINGMANIP_HPP

std::vector<std::string> split( const std::string &, const std::string & );

template <class T>
T str_to( const std::string &s) {
   /*Converts a string to numeric type T.*/

   T num;
   std::stringstream ss;
   
   ss << s;
   ss >> num;
   
   return num;
}

#endif
