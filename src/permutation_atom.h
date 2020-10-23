/*****************************************************************************

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

	    Header file for permutation_atom.cpp

******************************************************************************/

#ifndef __PERMUTATION_ATOM
#define __PERMUTATION_ATOM

#include <vector>
#include <set>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"
#include "include/label.hpp"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<int> > vector2_i;

class Permutation_atom{

    public: 

    Permutation_atom();
    ~Permutation_atom();

    imatrix permutation_lattice_translation_multilattice(const imatrix& snf, 
            const int& n_atom_primitive);

    imatrix permutation_lattice_translation(const imatrix& snf);

    imatrix permutation_rotation(
            const imatrix& snf, const imatrix& left, 
            const std::vector<imatrix>& rotate_matrix_array,
            const std::vector<dvector>& trans_vector_array,
            const std::vector<dvector>& position_primitive,
            const double& symprec);

    imatrix combine_table(const imatrix& rotation, 
        const imatrix& trans, const imatrix& snf);

    imatrix cluster_rotation(
            const std::vector<dvector>& cluster_candidate,
            const std::vector<dmatrix>& rotate_matrix_array, 
            const std::vector<dvector>& trans_vector_array,
            const std::vector<dvector>& position,
            const double& symprec,
            const dmatrix& inverse_axis_change);
};

#endif
