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

	    Class for converting derivative.out to structures
	
*****************************************************************************/

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "input.h"
#include "output.h"

#include "include/math.hpp"
#include "include/label.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<int> > vector2_i;

int main(){

    Input ip;
    ip.input_derivative_to_structure();
    dmatrix axis = ip.get_axis_primitive();
    std::vector<int> n_atoms = ip.get_n_atoms_primitive();
    std::vector<dvector> position_primitive = ip.get_position_primitive();
    std::vector<int> type = ip.get_type_primitive();
    const std::vector<int> n_type(ip.get_n_type());

    std::vector<imatrix> hnf_array = ip.get_hnf_array();
    std::vector<imatrix> snf_array = ip.get_snf_array();
    std::vector<imatrix> left_array = ip.get_left_array();
    vector2_i label_array = ip.get_label_array();
    int n_atom_label = label_array[0].size();

    std::cout.precision(10);
    std::cout.setf(std::ios::showpoint);

    //std::vector<dmatrix> new_axis_array;
//    vector2_i new_type_array;
 //   std::vector<std::vector<dvector> > new_position_array;
  //  std::vector<dmatrix> axis_change_array;
    std::vector<int> order;
    for (int i = 0; i < snf_array.size(); ++i){
        dmatrix new_axis;
        std::vector<dvector> new_position;
        std::vector<int> new_type;

        imatrix hnf = hnf_array[i];
        imatrix snf = snf_array[i];
        imatrix left = left_array[i];
        std::vector<int> label = label_array[i];

        dmatrix inverse_left;
        math::invert(left, inverse_left);

        // new_axis = axis * H * right = axis * (left)^-1 * S
        dmatrix axis_change = ublas::prod(inverse_left,snf);
        new_axis = ublas::prod(axis, axis_change);
 
        // lattice_label = S * z (z = fractional coodinate in new axis)
        std::vector<ivector> position_lattice;
        int n_lattice = snf(0,0) * snf(1,1) * snf(2,2);
        for (int j = 0; j < n_lattice; ++j){
            ivector lattice_label(3);
            int atom_index;
            sequential_number_to_lattice_expression
                (j, snf, lattice_label, atom_index);
            position_lattice.push_back(lattice_label);
        }

        // new_fractional_position = S^-1 * left * frac + S^-1 * lattice_label
        dmatrix inverse_snf;
        math::invert(snf, inverse_snf);
        int n(0);
        std::multimap<int, dvector> sort_position;
        for (int j = 0; j < position_primitive.size(); ++j){
            dvector internal_label = ublas::prod(left, position_primitive[j]);
            for (int k = 0; k < 3; ++k){
                while (internal_label(k) > 1 - 1e-5){
                    internal_label(k) -= 1.0;
                }
                while (internal_label(k) < -1e-5){
                    internal_label(k) += 1.0;
                }
            }
            for (int k = 0; k < n_lattice; ++k){
                dvector pos_in_newaxis = internal_label + position_lattice[k];
                dvector frac = ublas::prod(inverse_snf, pos_in_newaxis);
                int type_tmp;
                if (type[j] < n_type.size()){
                    type_tmp = label_array[i][n];
                }
                else {
                    type_tmp = type[j] + 1000;
                }
                std::pair<int, dvector> tmp_pair;
                tmp_pair.first = type_tmp;
                tmp_pair.second = frac;
                sort_position.insert(tmp_pair);
                ++n;
            }
        }
        std::multimap<int, dvector>::iterator it = sort_position.begin();
        while (it != sort_position.end()){
            new_type.push_back((*it).first);
            new_position.push_back((*it).second);
            ++it;
        }

        if (i == 0){
            std::set<int> tmp;
            for (int j = 0; j < new_type.size(); ++j){
                tmp.insert(new_type[j]);
            }
            order.resize(tmp.size());
            std::copy(tmp.begin(), tmp.end(), order.begin());
        }

        std::string str("structure"), str2, str3;
        std::stringstream ss, ss2;
        ss << n_atom_label;
        ss >> str2;
        str += str2;
        str += "_";
        ss2 << i + 1;
        ss2 >> str3;
        str += str3;

        std::ofstream output;
        output.open(str.c_str(), std::ios::out);

        Output_structure structure_out;
        structure_out.output_from_cell_to_file(output, new_axis, new_type,
                order, new_position);
        output << " axis_change = " << axis_change << std::endl;
        output << " n_type = ";
        for (int j = 0; j < n_type.size(); ++j){
            output << n_type[j] << " ";
        }
        output << std::endl;

//        new_axis_array.push_back(new_axis);
 //       new_type_array.push_back(new_type);
  //      new_position_array.push_back(new_position);
   //     axis_change_array.push_back(axis_change);
    }
  
    //Output op;
//    op.output_derivative_structure(new_axis_array, new_type_array, 
 //       new_position_array, axis_change_array, n_type);

    std::cout << " " << snf_array.size() 
        << " structures are generated." << std::endl;

    return 0;
}
