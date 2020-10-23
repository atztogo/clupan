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
	
    	Header file for sqs_binary.cpp

***************************************************************************/

#ifndef __SQS_BINARY
#define __SQS_BINARY

#include <iostream>
#include <vector>
#include <set>

#include <gsl/gsl_rng.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <correlation.h>

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<int> ivector;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;
typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<std::vector<int> > > vector3_i;
typedef std::vector<std::vector<std::vector<std::vector<int> > > > vector4_i;

class SQS_binary{

    std::vector<int> spin;
    dvector correlation_array, n_cluster_array;
    vector4_i coordination;
    vector3_i n_cluster_coordination;

    double disorder_cf;
    int n_step_ann;
    std::vector<double> temp_array;
    std::string criterion;

    void set_atomic_coordination(std::vector<imatrix>& permutation_array);
    double calc_spin_binary(const std::vector<int>& spin,
            const vector2_i& atom_cluster_coordination,
            const std::vector<int>& atom_n_cluster_coordination);

    public:

    SQS_binary(const std::vector<int>& spin_input,
            const dvector& correlation_array_input,      
            const dvector& n_cluster_array_input,
            std::vector<imatrix>& permutation_array_input,
            const int& n_step_ann_input, 
            const std::vector<double>& temp_array_input,
            const double& disorder_cf_input,
            const std::string& criterion);
    ~SQS_binary();

    void cmc();
    const dvector& get_correlation_array() const;
    const std::vector<int>& get_spin() const;
};

#endif
