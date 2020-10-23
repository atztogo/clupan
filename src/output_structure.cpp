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

	    Class for outputting structure
	
*****************************************************************************/

#include "output_structure.h"

Output_structure::Output_structure(){}

Output_structure::~Output_structure(){}

void Output_structure::output_from_structure(const dmatrix& axis, 
    const std::vector<int>& n_atoms, const std::vector<dvector>& position){

    dvector axis1(3), axis2(3), axis3(3);
    for (int i = 0; i < 3; ++i){
        axis1(i) = axis(i,0);
        axis2(i) = axis(i,1);
        axis3(i) = axis(i,2);
    }

    std::cout.precision(15);
    std::cout.setf(std::ios::showpoint);
    std::cout << " axis1 = " << axis1 << std::endl;
    std::cout << " axis2 = " << axis2 << std::endl;
    std::cout << " axis3 = " << axis3 << std::endl;

    std::cout << " n_atoms = "; 
    for (int i = 0; i < n_atoms.size(); ++i){
        std::cout << n_atoms[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < position.size(); ++i){
        std::cout << " position" << i+1 << " = " << position[i] << std::endl;
    }

}

void Output_structure::output_from_structure_to_file
    (std::ofstream& structure_out,
    const dmatrix& axis, const std::vector<int>& n_atoms, 
    const std::vector<dvector>& position){

    dvector axis1(3), axis2(3), axis3(3);
    for (int i = 0; i < 3; ++i){
        axis1(i) = axis(i,0);
        axis2(i) = axis(i,1);
        axis3(i) = axis(i,2);
    }

    structure_out.precision(15);
    structure_out.setf(std::ios::showpoint);
    structure_out << " axis1 = " << axis1 << std::endl;
    structure_out << " axis2 = " << axis2 << std::endl;
    structure_out << " axis3 = " << axis3 << std::endl;

    structure_out << " n_atoms = "; 
    for (int i = 0; i < n_atoms.size(); ++i){
        structure_out << n_atoms[i] << " ";
    }
    structure_out << std::endl;

    for (int i = 0; i < position.size(); ++i){
        structure_out << " position" << i+1 
            << " = " << position[i] << std::endl;
    }

}

void Output_structure::output_from_cell(const dmatrix& axis, 
    const std::vector<int>& type, 
    const std::vector<dvector>& position){

    dvector axis1(3), axis2(3), axis3(3);
    for (int i = 0; i < 3; ++i){
        axis1(i) = axis(i,0);
        axis2(i) = axis(i,1);
        axis3(i) = axis(i,2);
    }

    std::cout.precision(15);
    std::cout.setf(std::ios::showpoint);
    std::cout << " axis1 = " << axis1 << std::endl;
    std::cout << " axis2 = " << axis2 << std::endl;
    std::cout << " axis3 = " << axis3 << std::endl;

    std::vector<int> type_unique;
    for (int i = 0; i < type.size(); ++i){
        if (find(type_unique.begin(), type_unique.end(), type[i]) 
            == type_unique.end()){
            type_unique.push_back(type[i]);
        }
    }

    std::vector<int> n_atoms;
    int n(0);
    for (int i = 0; i < type_unique.size(); ++i){
        int n_i(0);
        for (int j = 0; j < type.size(); ++j){
            if (type_unique[i] == type[j]){
                std::cout << " position" << n+1 
                    << " = " << position[j] << std::endl;
                ++n;
                ++n_i;
            }
        }
        n_atoms.push_back(n_i);
    }

    std::cout << " n_atoms = "; 
    for (int i = 0; i < n_atoms.size(); ++i){
        std::cout << n_atoms[i] << " ";
    }
    std::cout << std::endl;

}

void Output_structure::output_from_cell_to_file(std::ofstream& structure_out,
    const dmatrix& axis, const std::vector<int>& type, 
    const std::vector<int>& order, const std::vector<dvector>& position){

    dvector axis1(3), axis2(3), axis3(3);
    for (int i = 0; i < 3; ++i){
        axis1(i) = axis(i,0);
        axis2(i) = axis(i,1);
        axis3(i) = axis(i,2);
    }

    structure_out.precision(15);
    structure_out.setf(std::ios::showpoint);
    structure_out << " axis1 = " << axis1 << std::endl;
    structure_out << " axis2 = " << axis2 << std::endl;
    structure_out << " axis3 = " << axis3 << std::endl;
/*`
    std::vector<int> type_unique;
    for (int i = 0; i < type.size(); ++i){
        if (find(type_unique.begin(), type_unique.end(), type[i]) 
            == type_unique.end()){
            type_unique.push_back(type[i]);
        }
    }
*/

    std::vector<int> n_atoms;
    int n(0);
    for (int i = 0; i < order.size(); ++i){
        int n_i(0);
        for (int j = 0; j < type.size(); ++j){
            if (order[i] == type[j]){
                structure_out << " position" << n+1 
                    << " = " << position[j] << std::endl;
                ++n;
                ++n_i;
            }
        }
        n_atoms.push_back(n_i);
    }

    structure_out << " n_atoms = "; 
    for (int i = 0; i < n_atoms.size(); ++i){
        structure_out << n_atoms[i] << " ";
    }
    structure_out << std::endl;

}
