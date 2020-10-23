/**************************************************************************

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
        
        Class for grand-canonical monte carlo simulation in binary system

***************************************************************************/

/*
    protected variables in class MC_common
    double k_b;
    gsl_rng * r;

    std::vector<int> spin, spin_array;
    dvector eci_ex_empty, correlation_array, n_cluster_array;

    int n_step_ann, n_step_eqv;
    std::vector<double> temp_array;
    std::vector<int> sublattice;

    double energy, energy_ave, energy_square, specific_heat;
    dvector correlation_ave;
*/

#include "gcmc_binary.h"

GCMC_binary::GCMC_binary(const Input& ip, const MC_initialize& mc_init)
    : MC_common(ip, mc_init){}

GCMC_binary::~GCMC_binary(){}

void GCMC_binary::gcmc(const Input& ip, const MC_initialize& mc_init){

    std::cout << "  energy of initial structure = " << energy << std::endl;

    // mu = mu(B) - mu(A)
    double mu = ip.get_mu()[1];

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, spin1, spin2, step_sum;
    double temperature, mu_change, diff_energy, prob, cut;
    dvector diff_correlation;

    energy_ave = 0.0; energy_square = 0.0;
    correlation_ave = ublas::zero_vector<double> (correlation_array.size());
    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_eqv, temperature, n_step);
        step_sum = 0;
        while (step_sum < n_step){
            std::vector<int> atoms;
            atoms.push_back(gsl_rng_uniform_int(r, n_atom_simulation));
            spin1 = spin[atoms[0]];
            mu_change = set_mu_change(spin_array, mu, spin1, spin2);
            diff_correlation 
                = gcmc_diff_correlation(spin, atoms, spin1, spin2, mc_init);
            diff_energy = ublas::inner_prod(eci_ex_empty, diff_correlation);
            prob = exp((mu_change - diff_energy)/ (k_b * temperature));
            cut =  gsl_rng_uniform(r);
            if (prob > cut){
                energy += diff_energy;
                correlation_array += diff_correlation;
                spin[atoms[0]] = spin2;
            }
            if (i == temp_array.size()){
                energy_ave += energy;
                energy_square += pow(energy, 2);
                correlation_ave += correlation_array;
            }
            ++step_sum;
        }
    }
    energy_ave /= double(n_step);
    energy_square /= double(n_step);
    correlation_ave /= double(n_step);
    specific_heat = (energy_square - pow (energy_ave,2)) 
        / k_b / pow(temperature, 2);
}

void GCMC_binary::gcmc_pseudo(const Input& ip, const MC_initialize& mc_init){

    std::cout << "  energy of initial structure = " << energy << std::endl;

    // mu = mu(B) - mu(A)
    double mu = ip.get_mu()[1];

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, spin1, spin2, step_sum;
    double temperature, mu_change, diff_energy, prob, cut;
    dvector diff_correlation;

    energy_ave = 0.0; energy_square = 0.0;
    correlation_ave = ublas::zero_vector<double> (correlation_array.size());
    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_eqv, temperature, n_step);
        step_sum = 0;
        while (step_sum < n_step){
            std::vector<int> atoms;
            atoms.push_back(gsl_rng_uniform_int(r, n_atom_simulation));
            spin1 = spin[atoms[0]];
            mu_change = set_mu_change(spin_array, mu, spin1, spin2);
            diff_correlation 
                = gcmc_diff_correlation(spin, atoms, spin1, spin2, mc_init);
            diff_energy = ublas::inner_prod(eci_ex_empty, diff_correlation);
            prob = exp((mu_change - diff_energy)/ (k_b * temperature));
            cut =  gsl_rng_uniform(r);
            if (prob > cut){
                energy += diff_energy;
                correlation_array += diff_correlation;
                spin[atoms[0]] = spin2;
            }
            if (i == temp_array.size()){
                energy_ave += energy;
                energy_square += pow(energy, 2);
                correlation_ave += correlation_array;
            }
            ++step_sum;
        }
    }
    energy_ave /= double(n_step);
    energy_square /= double(n_step);
    correlation_ave /= double(n_step);
    specific_heat = (energy_square - pow (energy_ave,2)) 
        / k_b / pow(temperature, 2);
}


dvector GCMC_binary::gcmc_diff_correlation
(std::vector<int>& spin, const std::vector<int>& atoms, 
 const int& spin1, const int& spin2,
 const MC_initialize& mc_init){

    dvector diff_correlation;

    dvector spin_old = calc_spin_binary
        (spin, atoms, mc_init.get_coordination(), 
         mc_init.get_n_cluster_coordination());
    spin[atoms[0]] = spin2;
    dvector spin_new = calc_spin_binary
        (spin, atoms, mc_init.get_coordination(), 
         mc_init.get_n_cluster_coordination());
    spin[atoms[0]] = spin1;

    diff_correlation = spin_new - spin_old;
    for (int i = 0; i < diff_correlation.size(); ++i){
        diff_correlation(i) /= mc_init.get_n_cluster_array()(i);
    }

    return diff_correlation;
}
