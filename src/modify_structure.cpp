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

	    Class for swapping atoms in structure file
	
*****************************************************************************/

#include <iostream>
#include <vector>

#include "getopt.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_input.h"
#include "parse_structure.h"
#include "output_structure.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc, char *argv[]){

    const char *file_name, *input_file_name = NULL;
    std::vector<int> sequence, combine, remove, empty;

    int option;
    opterr = 0;
    while((option=getopt(argc,argv,"f:i:"))!=-1){
        switch( option ){
            case 'f' : file_name = optarg;
                       break;
            case 'i' : input_file_name = optarg;
                       break;
        }
    }
    if (file_name == NULL){
        std::cerr << " Specify file name of structure (option -f)." 
            << std::endl;
        exit(8);
    }
    if (input_file_name == NULL){
        std::cerr << " Specify file name of input parameters (option -i)." 
            << std::endl;
        exit(8);
    }

    Parse_input in(input_file_name);
    in.assign("swap", sequence, empty);
    in.assign("combine", combine, empty);
    in.assign("remove", remove, empty);

    if (sequence.size() > 0 and combine.size() > 0){
        std::cerr << " Do not use 'swap' and 'combine' simultaneously." 
            << std::endl;
        exit(8);
    }
    if (sequence.size() > 0 and remove.size() > 0){
        std::cerr << " Do not use 'swap' and 'remove' simultaneously." 
            << std::endl;
        exit(8);
    }

    Parse_structure input(file_name);
    dmatrix axis(input.get_axis());
    std::vector<int> n_atoms(input.get_num_atoms());
    std::vector<dvector> position(input.get_position());

    std::vector<dvector> position_new;
    std::vector<int> n_atoms_new;

    if (sequence.size() > 0){
        for (int i = 0; i < sequence.size(); ++i){
            int index = sequence[i]-1;
            n_atoms_new.push_back(n_atoms[index]);
            int sum(0);
            for (int j = 0; j < index; ++j){
                sum += n_atoms[j];
            }
            for (int j = sum; j < sum + n_atoms[index]; ++j){
                position_new.push_back(position[j]);
            }
        }
        position = position_new;
        n_atoms = n_atoms_new;
    }

    if (remove.size() > 0){
        for (int i = 0; i < combine.size(); ++i){
            int new_index = combine[i];
            for (int j = 0; j < remove.size(); ++j){
                if (remove[j] < combine[i]){
                    --new_index;
                }
            }
            combine[i] = new_index;
        }
        int n(0);
        for (int i = 0; i < n_atoms.size(); ++i){
            if (find(remove.begin(), remove.end(), i+1) == remove.end()){
                n_atoms_new.push_back(n_atoms[i]);
                for (int j = n; j < n+n_atoms[i]; ++j){
                    position_new.push_back(position[j]);
                }
            }
            n += n_atoms[i];
        }
        position = position_new;
        n_atoms = n_atoms_new;
    }
    if (combine.size() > 0){
        int n_atoms1(0);
        std::vector<int> n_atoms2;
        std::vector<dvector> position1, position2;
        int n(0);
        for (int i = 0; i < n_atoms.size(); ++i){
            if (find(combine.begin(), combine.end(), i+1) != combine.end()){
                n_atoms1 += n_atoms[i];
                for (int j = n; j < n+n_atoms[i]; ++j){
                    position1.push_back(position[j]);
                }
            }
            else {
                n_atoms2.push_back(n_atoms[i]);
                for (int j = n; j < n+n_atoms[i]; ++j){
                    position2.push_back(position[j]);
                }
            }
            n += n_atoms[i];
        }
        position = position1;
        n_atoms.clear();
        n_atoms.push_back(n_atoms1);
        for (int i = 0; i < n_atoms2.size(); ++i){
            n_atoms.push_back(n_atoms2[i]);
        }
        for (int i = 0; i < position2.size(); ++i){
            position.push_back(position2[i]);
        }
    }
    
//    for (int i = 0; i < n_atoms_new.size(); ++i){
//        std::cout << n_atoms_new[i] << std::endl;
//    }
//    for (int i = 0; i < position_new.size(); ++i){
//        std::cout << position_new[i] << std::endl;
//    }

    Output_structure op;
    op.output_from_structure(axis, n_atoms, position);

    return 0;
}
