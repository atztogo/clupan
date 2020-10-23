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

	    Header file for sym.cpp

******************************************************************************/

#ifndef __SYM
#define __SYM

#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/math.hpp"

extern "C" int spg_get_symmetry(int rotation[][3][3], double translation[][3],
        int max_size, double lattice[3][3], double position[][3],
        int types[], int num_atom, double symprec);

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

class Sym{

    std::vector<imatrix> rotate_matrix_array;
    std::vector<dmatrix> double_rotate_matrix_array;
    std::vector<dvector> trans_vector_array;

    void get_sym(const dmatrix& axis, const std::vector<dvector>& position,
            const std::vector<int>& type, double symprec);

    public: 

    Sym(const dmatrix& axis, const std::vector<dvector>& position,
            const std::vector<int>& type, double symprec);

    ~Sym();

    void lattice_based_primitive(const dmatrix& axis, 
            const std::vector<dvector>& position, 
            const std::vector<int>& type_lattice, const double& symprec,
            const std::vector<int> n_atoms, const dmatrix& axis_change, 
            const dmatrix& inverse_axis_change);

    void supercell_symmetry(const dmatrix& axis_change);

    const std::vector<imatrix>& get_rotate_matrix() const;
    const std::vector<dvector>& get_trans_vector() const;
    const std::vector<dmatrix>& get_double_rotate_matrix() const;

};

#endif
