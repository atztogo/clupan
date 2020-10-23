/****************************************************************************

        Copyright (C) 2013 Atsuto Seko
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

	    Header file for mc_common.cpp
		
****************************************************************************/

#ifndef __MC_COMMON
#define __MC_COMMON

#include <vector>
#include <set>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "gsl/gsl_rng.h"

#include "input.h"
#include "mc_initialize.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<int> > vector2_i;
typedef std::vector<std::vector<std::vector<int> > > vector3_i;
typedef std::vector<std::vector<std::vector<std::vector<int> > > > vector4_i;

class MC_common{

    protected:

    double k_b;
    gsl_rng * r;

    std::vector<int> spin, spin_array;
    dvector eci_ex_empty, correlation_array, n_cluster_array;

    int n_step_ann, n_step_eqv;
    std::vector<double> temp_array;
    std::vector<int> sublattice;

    double energy, energy_ave, energy_square, specific_heat;
    dvector correlation_ave;

    double get_sqs_score
        (const dvector& correlation_array_new,
         const dvector& disorder_cf, const std::string& criterion);

    double set_mu_change
        (const std::vector<int>& spin_array, const double& mu, 
         const int& spin1, int& spin2);
    double set_mu_change
        (const std::vector<int>& spin_array, const std::vector<double>& mu, 
         const int& spin1, int& spin2);

    void set_temp_and_n_step
        (const int& temp_step, const std::vector<double>& temp_array, 
         const int& n_step_ann, const int& n_step_eqv, 
         double& temperature, int& n_step);

    void swap_spin(std::vector<int>& spin, const int& atom1, const int& atom2);
    void swap_charge
        (std::vector<double>& charge, const int& atom1, const int& atom2);

    dvector calc_spin_binary
        (const std::vector<int>& spin, const std::vector<int>& atoms,
         const vector4_i& coordination, 
         const vector3_i& n_cluster_coordination);

    dvector calc_spin_cmc
        (const std::vector<int>& spin, const int& atom1, const int& atom2,
         const vector4_i& coordination, 
         const std::vector<dmatrix>& cluster_function_array);

    dvector calc_spin_gcmc
        (const std::vector<int>& spin, const int& atom1,
         const vector4_i& coordination, 
         const std::vector<dmatrix>& cluster_function_array);


    private:

    double calc_local_spin_binary
        (const std::vector<int>& spin,
         const vector2_i& atom_cluster_coordination,
         const std::vector<int>& atom_n_cluster_coordination);

    double calc_local_spin
        (const std::vector<int>& spin,
         const vector2_i& atom_cluster_coordination,
         const dmatrix& cluster_function);

    double calc_local_spin
        (const std::vector<int>& spin, const int& other_atom,
         const vector2_i& atom_cluster_coordination,
         const dmatrix& cluster_function);

    double calc_product
        (const std::vector<int>& spin_values, const dmatrix& function);

    public:

    MC_common(const Input& ip, const MC_initialize& mc_init);
    ~MC_common();


    const double& get_energy() const;
    const dvector& get_correlation_array() const;
    const std::vector<int>& get_spin() const;

    const double& get_energy_average() const;
    const double& get_energy_square() const;
    const dvector& get_correlation_average() const;
    const double& get_specific_heat() const;

};

#endif
