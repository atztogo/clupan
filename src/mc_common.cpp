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

	    Common functions for MC simulations
		
****************************************************************************/

#include "mc_common.h"

MC_common::MC_common(const Input& ip, const MC_initialize& mc_init){

    k_b = 0.00008617343;

    unsigned long init = time(NULL) + getpid();
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, init);

    spin_array = ip.get_spin_array();
    temp_array = ip.get_temp_array();

    spin = mc_init.get_spin();
    n_step_ann = mc_init.get_n_step_ann();
    n_step_eqv = mc_init.get_n_step_eqv();

    eci_ex_empty = mc_init.get_eci_without_empty();
    energy = mc_init.get_energy();
    correlation_array = mc_init.get_correlation_array();
    n_cluster_array = mc_init.get_n_cluster_array();
    sublattice = mc_init.get_sublattice();
}

MC_common::~MC_common(){}

double MC_common::get_sqs_score
(const dvector& correlation_array_new,
const dvector& disorder_cf, const std::string& criterion){

    double score_new(0.0);
    int n_cluster = correlation_array_new.size();
    dvector diff = correlation_array_new - disorder_cf;
    if (criterion == "rms"){
        score_new = ublas::inner_prod(diff, diff);
    }
    else if (criterion == "abs"){
        for (int i = 0; i < n_cluster; ++i){
            score_new += fabs(diff(i));
        }
    }

    return score_new;
}
// available only for single lattice CE
double MC_common::set_mu_change
(const std::vector<int>& spin_array, const double& mu, 
 const int& spin1, int& spin2){

    double mu_change;
    if (spin1 == spin_array[0]){
        spin2 = spin_array[1];
        mu_change = mu;

    }
    else if (spin1 == spin_array[1]){
        spin2 = spin_array[0];
        mu_change = -mu;
    }

    return mu_change;
}

// available only for single lattice CE
double MC_common::set_mu_change
(const std::vector<int>& spin_array, const std::vector<double>& mu, 
 const int& spin1, int& spin2){

    double mu_change;
    std::vector<int> candidate;
    int orig; 
    for (int i = 0; i < spin_array.size(); ++i){
        if (spin1 != spin_array[i]){
            candidate.push_back(i);
        }
        else {
            orig = i;
        }
    }
    int index_new = candidate[rand()%candidate.size()];
    spin2 = spin_array[index_new];

    mu_change = mu[index_new] - mu[orig];

    return mu_change;
}

void MC_common::set_temp_and_n_step
(const int& temp_step, const std::vector<double>& temp_array, 
 const int& n_step_ann, const int& n_step_eqv, 
 double& temperature, int& n_step){

    if (temp_step < temp_array.size()){
        temperature = temp_array[temp_step];
        n_step = n_step_ann;
    }
    else {
        temperature = temp_array[temp_array.size()-1];
        n_step = n_step_eqv;
    }
    std::cout << "  temperature = " << temperature << std::endl;
    std::cout << "  n_step = " << n_step << std::endl;
}

void MC_common::swap_spin
(std::vector<int>& spin, const int& atom1, const int& atom2){

    std::vector<int>::iterator it1, it2;
    it1 = spin.begin() + atom1;
    it2 = spin.begin() + atom2;
    iter_swap(it1, it2);
}

void MC_common::swap_charge
(std::vector<double>& charge, const int& atom1, const int& atom2){

    std::vector<double>::iterator it1, it2;
    it1 = charge.begin() + atom1;
    it2 = charge.begin() + atom2;
    iter_swap(it1, it2);
}

dvector MC_common::calc_spin_binary
(const std::vector<int>& spin, const std::vector<int>& atoms,
 const vector4_i& coordination, const vector3_i& n_cluster_coordination){

    int n_cluster = coordination[atoms[0]].size();
    dvector spin_around_atom = ublas::zero_vector<double> (n_cluster);
    for (int i = 0; i < atoms.size(); ++i){
        for (int j = 0; j < n_cluster; ++j){
            spin_around_atom(j) += calc_local_spin_binary
                (spin, coordination[atoms[i]][j], 
                 n_cluster_coordination[atoms[i]][j]);
        }
    }

    return spin_around_atom;
}
// local function
double MC_common::calc_local_spin_binary
(const std::vector<int>& spin,
 const vector2_i& atom_cluster_coordination,
 const std::vector<int>& atom_n_cluster_coordination){

    double spin_sum(0.0);
    for (int i = 0; i < atom_cluster_coordination.size(); ++i){
        double product(1.0);
        for (int j = 0; j < atom_cluster_coordination[i].size(); ++j){
            int atom = atom_cluster_coordination[i][j];
            product *= double(spin[atom]);
        }
        spin_sum += product * atom_n_cluster_coordination[i];
    }

    return spin_sum;
}

dvector MC_common::calc_spin_cmc
(const std::vector<int>& spin, const int& atom1, const int& atom2,
 const vector4_i& coordination, 
 const std::vector<dmatrix>& cluster_function_array){

    int n_cluster = coordination[atom1].size();
    dvector spin_around_atom = ublas::zero_vector<double> (n_cluster);

    for (int i = 0; i < n_cluster; ++i){
        spin_around_atom(i) += calc_local_spin
            (spin, coordination[atom1][i], cluster_function_array[i]);
    }
    for (int i = 0; i < n_cluster; ++i){
        spin_around_atom(i) += calc_local_spin
            (spin, atom1, coordination[atom2][i], cluster_function_array[i]);
    }

    return spin_around_atom;
}

dvector MC_common::calc_spin_gcmc
(const std::vector<int>& spin, const int& atom1,
 const vector4_i& coordination, 
 const std::vector<dmatrix>& cluster_function_array){

    int n_cluster = coordination[atom1].size();
    dvector spin_around_atom = ublas::zero_vector<double> (n_cluster);

    for (int i = 0; i < n_cluster; ++i){
        spin_around_atom(i) += calc_local_spin
            (spin, coordination[atom1][i], cluster_function_array[i]);
    }

    return spin_around_atom;
}


// local function
double MC_common::calc_local_spin
(const std::vector<int>& spin,
 const vector2_i& atom_cluster_coordination,
 const dmatrix& cluster_function){

    double spin_sum(0.0);
    for (int i = 0; i < atom_cluster_coordination.size(); ++i){
        std::vector<int> spin_values;
        for (int j = 0; j < atom_cluster_coordination[i].size(); ++j){
            int atom = atom_cluster_coordination[i][j];
            spin_values.push_back(spin[atom]);
        }
        double product = calc_product(spin_values, cluster_function);
        spin_sum += product;
    }

    return spin_sum;
}

// local function
double MC_common::calc_local_spin
(const std::vector<int>& spin, const int& other_atom,
 const vector2_i& atom_cluster_coordination,
 const dmatrix& cluster_function){

    double spin_sum(0.0);
    for (int i = 0; i < atom_cluster_coordination.size(); ++i){
        if (find(atom_cluster_coordination[i].begin(), 
                    atom_cluster_coordination[i].end(), other_atom) 
                == atom_cluster_coordination[i].end()){
            std::vector<int> spin_values;
            for (int j = 0; j < atom_cluster_coordination[i].size(); ++j){
                int atom = atom_cluster_coordination[i][j];
                spin_values.push_back(spin[atom]);
            }
            double product = calc_product(spin_values, cluster_function);
            spin_sum += product;
        }
    }

    return spin_sum;
}

// local function
double MC_common::calc_product
(const std::vector<int>& spin_values, const dmatrix& function){

    double product(1.0);
    for (int i = 0; i < function.size1(); ++i){
        double sum(0.0);
        for (int j = 0; j < function.size2(); ++j){
            double tmp_prod(1.0);
            for (int k = 0; k < j; ++k){
                tmp_prod *= spin_values[i];
            }
            sum += function(i,j) * tmp_prod;
        }
        product *= sum;
    }

    return product;
}

const double& MC_common::get_energy() const{
    return energy;
}
const std::vector<int>& MC_common::get_spin() const{
    return spin;
}
const dvector& MC_common::get_correlation_array() const{
    return correlation_array;
}
const double& MC_common::get_energy_average() const{
    return energy_ave;
}
const dvector& MC_common::get_correlation_average() const{
    return correlation_ave;
}
const double& MC_common::get_energy_square() const{
    return energy_square;
}
const double& MC_common::get_specific_heat() const{
    return specific_heat;
}
