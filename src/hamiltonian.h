/******************************************************************************

        Copyright (C) 2007 Atsuto Seko
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


	Header file for hamiltonian.cpp

******************************************************************************/

#ifndef __HAMILTONIAN
#define __HAMILTONIAN

#include <iostream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

class Hamiltonian{

    dmatrix onsiteFunction, clusterFunction;

    dmatrix gramSchmidt(const std::vector<int>& spin);

    dvector functionProduct(const dvector& func1, 
            const dvector& func2);

    double innerProduct(const std::vector<int>& spin, 
            const dvector& func1, const dvector& func2);

    public: 

    Hamiltonian(const std::vector<int>& spin);
    ~Hamiltonian();

    const dmatrix& get_onsite_function() const;

};

#endif
