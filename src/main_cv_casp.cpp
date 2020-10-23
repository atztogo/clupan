/****************************************************************************

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


	Main program for least square fitting 

*****************************************************************************/

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "input.h"
#include "output.h"
#include "least_squares.h"
#include "cv.h"

namespace ublas = boost::numeric::ublas;
typedef boost::numeric::ublas::vector<double> dvector;
typedef boost::numeric::ublas::matrix<double> dmatrix;

int main(int argc, char *argv[]){

    Input ip;
    ip.input_cv_casp(argc, argv);
    const dvector energy(ip.get_energy());
    const dmatrix correlation(ip.get_correlation());
    const dmatrix weight(ip.get_weight());
    const std::vector<std::string> 
        energy_structure_index(ip.get_energy_structure_index());

    const std::vector<int> group_index(ip.get_group_index_array());
    const std::vector<std::string> 
        structure_index_array(ip.get_structure_index_array());

    int n_structure = correlation.size1();
    int n_cluster = correlation.size2();

    std::vector<int> group;
    for (int i = 0; i < energy.size(); ++i){
        for (int j = 0; j < structure_index_array.size(); ++j){
            if (structure_index_array[j] == energy_structure_index[i]){
                group.push_back(group_index[j]);
                break;
            }
            else if (j == structure_index_array.size() - 1){
                std::cerr << "  structure " << energy_structure_index[i] 
                    << " is not found in group file." << std::endl;
                exit(8);
            }
        }
    }

    double cv_score;
    std::vector<double> cv_casp, cv_casp2;
    std::vector<int> n_structure_in_group;

    Least_squares least;
    int check = least.fit(energy, correlation, weight);
    if (check == 1){
        dvector eci = least.get_eci();
        CV cv;
        cv_score = cv.calc_cv_score(energy, correlation, weight, eci);
        cv_casp = cv.calc_cv_casp(energy, correlation, eci, group);
        n_structure_in_group = cv.get_n_structure_in_group();
        Output op;
        op.output_cv_casp(cv_casp, n_structure_in_group, cv_score);
    }
    else {
        std::cout << " matrix of correlation functions is singular. " 
            << std::endl;
    }
	return 0;

}

