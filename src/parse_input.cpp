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

#include "parse_input.h"

Parse_input::Parse_input(const char *filename){

    std::ifstream input(filename);

    if (input.fail()){
        std::cerr << "Error: Could not open " << filename << "\n";
        exit(8);
    }

    std::vector<std::string> dataLine;

    std::string line;
    while (getline(input,line) && !input.eof()){
        int place = line.find("=");
        if (place > 0){
            line.replace(place, 1, " ");
            std::stringstream str;
            str << line;
            std::istream_iterator<std::string> iss_iter (str);
            std::istream_iterator<std::string> iss_end;
            for (; iss_iter != iss_end; ++iss_iter){
                dataLine.push_back(*iss_iter);
            }
            data.push_back(dataLine);
        }
        dataLine.clear();
    }
}

Parse_input::~Parse_input(){}

