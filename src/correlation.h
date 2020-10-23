/***************************************************************************

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
	
    	Header file for correlation.cpp

***************************************************************************/

#ifndef __CORRELATION
#define __CORRELATION

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "permutation.h"
#include "permutation_atom.h"
#include "hamiltonian.h"

typedef std::vector<std::vector<std::pair<ivector, int> > > vector2_pair;

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<int> ivector;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

class Correlation{

    double calc_product(std::vector<int>& spin_values, dmatrix& function);

    public:

    Correlation();
    ~Correlation();

    std::vector<imatrix> get_permutation_array(
            const vector2_pair& unique_cluster_array,
            const std::vector<dvector> position_primitive,
            const std::vector<dvector> position,
            const dmatrix& inverse_axis_change,
            const std::vector<dmatrix>& rotate_matrix_array,
            const std::vector<dvector>& trans_vector_array,
            const double& symprec);

    std::vector<dmatrix> get_cluster_function_array(
            const std::vector<int>& spin_array,
            std::vector<imatrix>& permutation_array,
            std::vector<int>& cluster_number_array);

    std::vector<dmatrix> get_cluster_function_array_multilattice(
            const std::vector<int>& n_type,
            const std::vector<int>& spin_array,
            const vector2_i& cluster_sublattice,
            std::vector<imatrix>& permutation_array,
            std::vector<int>& cluster_number_array);


    dvector get_correlation_function(
            const std::vector<imatrix>& permutation_array,
            const std::vector<dmatrix>& cluster_function_array,
            const std::vector<int>& spin);

};

#endif
