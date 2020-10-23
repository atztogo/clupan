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

#undef NDEBUG

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

int main(){

    Input ip;
    ip.input_lsf();
    std::vector<int> cluster_index = ip.get_cluster_index();
    dvector energy = ip.get_energy();
    dmatrix correlation = ip.get_correlation();
    dmatrix weight = ip.get_weight();

/*    dmatrix correlation_calc(n_structure, n_cluster);
    for (int i = 0; i < n_structure; ++i){
        for (int j = 0; j < n_cluster; ++j){
            correlation_calc(i,j) = correlation(i,cluster_index[j]);
        }
    }
*/
    Least_squares least;
    int check = least.fit(energy, correlation, weight);
    if (check == 1){
        dvector eci = least.get_eci();
        double square_error = least.get_square_error();
        CV cv;
        double cv_score = cv.calc_cv_score(energy, correlation, weight, eci);

        dmatrix prod_mat = ublas::prod(weight, correlation);
        dmatrix trans_mat = ublas::trans(correlation);
        dmatrix prod_mat2 = ublas::prod(trans_mat, prod_mat);
        dmatrix inverse;
        math::invert(prod_mat2, inverse);

        Output op;
        op.output_lsf(cluster_index, eci, square_error, 
            cv_score, inverse, check);
    }
    else {
        std::cout << " matrix of correlation functions is singular. " 
            << std::endl;
    }

	return 0;

}

