/**************************************************************************

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
        
        Class for monte carlo simulation in N-nary system

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

#include "mc.h"

MC::MC(const Input& ip, const MC_initialize& mc_init)
    : MC_common(ip, mc_init){}

MC::~MC(){}

void MC::find_sqs(const Input& ip, const MC_initialize& mc_init){

    const dvector disorder_cf = mc_init.get_disorder_cf();
    const std::string criterion = ip.get_criterion();

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, atom1, atom2, step_sum;
    double temperature, prob, cut, score, score_new, score_min, diff_score;
    dvector diff_correlation, correlation_array_min;
    std::vector<int> spin_min;

    score = get_sqs_score(correlation_array, disorder_cf, criterion);
    score_min = score;
    spin_min = spin;
    correlation_array_min = correlation_array;

    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_ann, temperature, n_step);
        step_sum = 0; 
        while (step_sum < n_step){
            atom1 = gsl_rng_uniform_int(r, n_atom_simulation);
            atom2 = gsl_rng_uniform_int(r, n_atom_simulation);
            if (spin[atom1] != spin[atom2] 
                    and sublattice[atom1] == sublattice[atom2]){
                diff_correlation 
                    = cmc_diff_correlation(spin, atom1, atom2, mc_init);
                score_new = get_sqs_score
                    (correlation_array + diff_correlation, 
                    disorder_cf, criterion);
                diff_score = score_new - score;
                prob = exp(- diff_score / (k_b * temperature));
                cut =  gsl_rng_uniform(r);
                if (prob > cut){
                    correlation_array += diff_correlation;
                    score = score_new;
                    swap_spin(spin, atom1, atom2);
                    if (score < score_min){
                        score_min = score;
                        spin_min = spin;
                        correlation_array_min = correlation_array;
                    }
                }
                ++step_sum;
            }
        }
        std::cout << "   rms = " 
            << sqrt(score/double(n_cluster)) << std::endl;
    }
    spin = spin_min;
    correlation_array = correlation_array_min;

    if (criterion == "rms"){
        std::cout << "   rms (min, sqs) = "
            << sqrt(score_min/double(n_cluster)) << std::endl;
    }
    else if (criterion == "abs"){
        std::cout << "   mean of abs (min, sqs) = "
            << score_min/double(n_cluster) << std::endl;
    }
}

void MC::gcmc(const Input& ip, const MC_initialize& mc_init){

    // mu = mu(B) - mu(A)
    std::vector<double> mu = ip.get_mu();

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, atom1, spin1, spin2, step_sum;
    double temperature, mu_change, diff_energy, prob, cut;
    dvector diff_correlation;

    energy_ave = 0.0; energy_square = 0.0;
    correlation_ave = ublas::zero_vector<double> (correlation_array.size());
    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_eqv, temperature, n_step);
        step_sum = 0;
        while (step_sum < n_step){
            atom1 = gsl_rng_uniform_int(r, n_atom_simulation);
            spin1 = spin[atom1];
            mu_change = set_mu_change(spin_array, mu, spin1, spin2);
            diff_correlation 
                = gcmc_diff_correlation(spin, atom1, spin2, mc_init);
            diff_energy = ublas::inner_prod(eci_ex_empty, diff_correlation);
            prob = exp((mu_change - diff_energy)/ (k_b * temperature));
            cut =  gsl_rng_uniform(r);
            if (prob > cut){
                energy += diff_energy;
                correlation_array += diff_correlation;
                spin[atom1] = spin2;
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

void MC::cmc(const MC_initialize& mc_init){

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, atom1, atom2, step_sum;
    double temperature, diff_energy, prob, cut;
    dvector diff_correlation;

    energy_ave = 0.0; energy_square = 0.0;
    correlation_ave = ublas::zero_vector<double> (correlation_array.size());
    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_eqv, temperature, n_step);
        step_sum = 0; 
        while (step_sum < n_step){
            atom1 = gsl_rng_uniform_int(r, n_atom_simulation);
            atom2 = gsl_rng_uniform_int(r, n_atom_simulation);
            if (spin[atom1] != spin[atom2] 
                    and sublattice[atom1] == sublattice[atom2]){
                diff_correlation 
                    = cmc_diff_correlation(spin, atom1, atom2, mc_init);
                diff_energy = ublas::inner_prod(eci_ex_empty, diff_correlation);
                prob = exp(- diff_energy / (k_b * temperature));
                cut =  gsl_rng_uniform(r);
                if (prob > cut){
                    energy += diff_energy;
                    correlation_array += diff_correlation;
                    swap_spin(spin, atom1, atom2);
                }
                if (i == temp_array.size()){
                    energy_ave += energy;
                    energy_square += pow(energy, 2);
                    correlation_ave += correlation_array;
                }
                ++step_sum;
            }
        }
    }
    energy_ave /= double(n_step);
    energy_square /= double(n_step);
    correlation_ave /= double(n_step);
    specific_heat = (energy_square - pow (energy_ave,2)) 
        / k_b / pow(temperature, 2);
}

void MC::cmc_ewald(const Input& ip, const MC_initialize& mc_init){

    const double epsilon = ip.get_epsilon();
    std::vector<double> charge = mc_init.get_charge();

    Ewald ewald_obj = mc_init.get_ewald_obj();
    double ewald_energy = ewald_obj.get_energy() / epsilon;
    std::cout << "  electrostatic energy of initial structure = ";
    std::cout.precision(10);
    std::cout << ewald_energy << std::endl;

    std::vector<std::complex<double> > charge_reciprocal
        = ewald_obj.initial_calc_charge_reciprocal(charge);

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    int n_step, atom1, atom2, step_sum;
    double temperature, diff_ewald, diff_energy, prob, cut;
    dvector diff_correlation;

    energy_ave = 0.0; energy_square = 0.0;
    correlation_ave = ublas::zero_vector<double> (correlation_array.size());
    for (int i = 0; i < temp_array.size() + 1; ++i){
        set_temp_and_n_step
            (i, temp_array, n_step_ann, n_step_eqv, temperature, n_step);
        step_sum = 0; 
        while (step_sum < n_step){
            atom1 = gsl_rng_uniform_int(r, n_atom_simulation);
            atom2 = gsl_rng_uniform_int(r, n_atom_simulation);
            if (spin[atom1] != spin[atom2] 
                    and sublattice[atom1] == sublattice[atom2]){
                diff_correlation 
                    = cmc_diff_correlation(spin, atom1, atom2, mc_init);
                diff_ewald = ewald_obj.calc_mc_energy
                    (charge, charge_reciprocal, atom1, atom2);
                diff_energy = ublas::inner_prod(eci_ex_empty, diff_correlation)
                    + diff_ewald / epsilon;
                prob = exp(- diff_energy / (k_b * temperature));
                cut =  gsl_rng_uniform(r);
                if (prob > cut){
                    energy += diff_energy;
                    ewald_energy += diff_ewald / epsilon;
                    correlation_array += diff_correlation;
                    ewald_obj.calc_diff_charge_reciprocal
                        (charge, charge_reciprocal, atom1, atom2);
                    swap_spin(spin, atom1, atom2);
                    swap_charge(charge, atom1, atom2);
                }
                if (i == temp_array.size()){
                    energy_ave += energy;
                    energy_square += energy * energy;
                    correlation_ave += correlation_array;
                }
                ++step_sum;
            }
        }
    }
    energy_ave /= double(n_step);
    energy_square /= double(n_step);
    correlation_ave /= double(n_step);
    specific_heat = (energy_square - pow (energy_ave,2)) 
        / k_b / pow(temperature, 2);
    std::cout << "  Electrostatic energy of final structure = "
        << ewald_energy << std::endl;
}

dvector MC::cmc_diff_correlation
(std::vector<int>& spin, const int& atom1, const int& atom2, 
 const MC_initialize& mc_init){

    dvector diff_correlation, spin_old, spin_new;

    spin_old = calc_spin_cmc
        (spin, atom1, atom2, mc_init.get_coordination(), 
         mc_init.get_cluster_function_array());
    swap_spin(spin, atom1, atom2);
    spin_new = calc_spin_cmc
        (spin, atom1, atom2, mc_init.get_coordination(), 
         mc_init.get_cluster_function_array());
    swap_spin(spin, atom1, atom2);

    diff_correlation = spin_new - spin_old;
    for (int i = 0; i < diff_correlation.size(); ++i){
        diff_correlation(i) /= mc_init.get_n_cluster_array()(i);
    }

    return diff_correlation;
}

dvector MC::gcmc_diff_correlation
(std::vector<int>& spin, const int& atom1, const int& spin2, 
 const MC_initialize& mc_init){

    dvector diff_correlation, spin_old, spin_new;

    int spin1 = spin[atom1];
    spin_old = calc_spin_gcmc
        (spin, atom1, mc_init.get_coordination(), 
         mc_init.get_cluster_function_array());
    spin[atom1] = spin2;
    spin_new = calc_spin_gcmc
        (spin, atom1, mc_init.get_coordination(), 
         mc_init.get_cluster_function_array());
    spin[atom1] = spin1;

    diff_correlation = spin_new - spin_old;
    for (int j = 0; j < diff_correlation.size(); ++j){
        diff_correlation(j) /= mc_init.get_n_cluster_array()(j);
    }

    return diff_correlation;
}
