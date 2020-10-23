/******************************************************************************

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

	    Class for parsing POSCAR
	
*****************************************************************************/

#include "parse_poscar.h"

Parse_poscar::Parse_poscar(const char *filename){

    std::ifstream input(filename);

    if (input.fail()){
        std::cerr << "Error: Could not open " << filename << "\n";
        exit(8);
    }	

    // line 1 : comment
    std::getline( input, comment );

    // local variables
    std::stringstream str;
    std::string line;

    // line 2 : length unit
    std::getline( input, line );
    str << line;
    double unit;
    str >> unit;
    str.clear();

    // line 3-5 : axis
    double x, y, z;
    axis = ublas::identity_matrix<double>(3);
    for (int i = 0; i < 3; ++i){
        std::stringstream str2;
        std::getline( input, line );
        str2 << line;
        str2 >> x >> y >> z;
        axis(0,i) = x * unit;
        axis(1,i) = y * unit;
        axis(2,i) = z * unit;
    }

    // line 6 : number of atoms 
    std::getline( input, line );
    std::getline( input, line );
    
    str << line;
    for (int i = 0; i < 20; ++i){
        int a = -1;
        str >> a;
        if ( a == -1)
            break;
        else if ( a != -1)
            num_atoms.push_back(a);
    }	
    str.clear();

    // line 7 : coordination type
    std::getline (input, coordinateType);

    // line 8- : position of atoms
    int num_atoms_all = 0;
    for (int i = 0; i < num_atoms.size(); ++i)
        num_atoms_all += num_atoms[i];

    std::string posall;
    for (int i = 0; i < num_atoms_all; ++i){
        std::stringstream str2;
        std::getline(input, posall);
        str2 << posall;
        str2 >> x >> y >> z;
        dvector tmp_position(3);
        tmp_position(0) = modify_position(x);
        tmp_position(1) = modify_position(y);
        tmp_position(2) = modify_position(z);
        position.push_back(tmp_position);
    }

}

Parse_poscar::~Parse_poscar(){}

const dmatrix& Parse_poscar::get_axis() const{

    return axis;
}

const std::vector<int>& Parse_poscar::get_num_atoms() const{

    return num_atoms;
}

const std::vector<dvector>& Parse_poscar::get_position() const{

    return position;
}

double Parse_poscar::modify_position(double x) const{

    double x_mod;

    if (x > 1 - 1e10)
        x_mod = x - floor(x);
    else if (x < - 1e10)
        x_mod = 1 + x - floor(x);
    else
        x_mod = x;

    return x_mod;

}


