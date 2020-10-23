/************************************************************************

        Copyright (C) 2011 Atsuto Seko
                seko@cms.mtl.kyoto-u.ac.jp

        This program is free software; you can redistribute it and/or
        modify it under the terms of the GNU General Public License
        as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to
        the Free Software Foundation, Inc., 51 Franklin Street,
        Fifth Floor, Boston, MA 02110-1301, USA, or see
        http://www.gnu.org/copyleft/gpl.txt
	

	    Header file for parsing input files
		
************************************************************************/

#ifndef __PARSEINPUT
#define __PARSEINPUT

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <iterator>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

class Parse_input{


    std::vector<std::vector<std::string> > data;

    public:
    Parse_input(const char *inputFile);
    ~Parse_input();

    template<typename T> void assign
        (const char* name, T& variable, T defaultValue); 
    template<typename T> void assign
        (const char* name, std::vector<T>& variable,
         std::vector<T>& defaultValue);
    template<typename T> void assignNeed
        (const char* name, T& variable);
    template<typename T> void assignNeed
        (const char* name, std::vector<T>& variable);

};

template<typename T> 
void Parse_input::assign(const char* name, T& variable, T defaultValue){

    const int data_size = data.size();
    int line(data_size);
    for (int i = 0; i < data_size; ++i){
        if (data[i][0] == name)
            line = i;
    }

    if (line < data.size())
        variable = boost::lexical_cast<T>(data[line][1]);
    else {
        variable = defaultValue;
        std::clog << "  default parameter : " << name << " = "
            << variable << "\n";
    }
}

template<typename T> 
void Parse_input::assign(const char* name, std::vector<T>& variable,
        std::vector<T>& defaultValue){

    const int data_size = data.size();
    int line(data_size);
    for (int i = 0; i < data_size; ++i){
        if (data[i][0] == name)
            line = i;
    }
    if (line < data.size()){
        for (int i = 1; i < data[line].size(); ++i){
            variable.push_back(boost::lexical_cast<T>(data[line][i]));
        }
    }
    else {
        variable = defaultValue;
        std::clog << "  default parameter : " << name << " = ";
        for (int i = 0; i < variable.size(); ++i){
            std::clog << variable[i] << " ";
        }
        std::clog << "\n";
    }
}

template<typename T> 
void Parse_input::assignNeed(const char* name, T& variable){

    const int data_size = data.size();
    int line(data_size);
    for (int i = 0; i < data_size; ++i){
        if (data[i][0] == name)
            line = i;
    }
    if (line < data.size())
        variable = boost::lexical_cast<T>(data[line][1]);
    else {
        std::cerr << "  Set the variable " << name << " !!" << std::endl;
        exit(8);
    }
}

template<typename T> 
void Parse_input::assignNeed(const char* name, std::vector<T>& variable){

    const int data_size = data.size();
    int line(data_size);
    for (int i = 0; i < data_size; ++i){
        if (data[i][0] == name)
            line = i;
    }
    if (line < data.size()){
        for (int i = 1; i < data[line].size(); ++i){
            variable.push_back(boost::lexical_cast<T>(data[line][i]));
        }
    }
    else {
        std::cerr << "  Set the variable " << name << " !!" << std::endl;
        exit(8);
    }
}

#endif
