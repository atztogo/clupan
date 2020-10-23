/******************************************************************************

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

	    Class for converting correlation for casp
	
*****************************************************************************/

#include <iostream>
#include <vector>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parse_input.h"
#include "input.h"
#include "ewald.h"
#include "least_squares.h"
#include "cv.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<double> dmatrix;
typedef std::vector<std::vector<double> > vector2_d;

vector2_d get_charge_permutation(const vector2_d& charge_range){

    vector2_d permutation;
    std::vector<double> candidate;
    std::vector<int> n_charge_array;
    for (int i = 0; i < charge_range.size(); ++i){
        for (int j = 0; j < charge_range[i].size(); ++j){
            candidate.push_back(charge_range[i][j]);
        }
        n_charge_array.push_back(charge_range[i].size());
    }

    long n_comb(1);
    std::vector<int> zero;
    for (int i = 0; i < charge_range.size(); ++i){
        n_comb *= charge_range[i].size();
        zero.push_back(0);
    }

    for (long i = 0; i < n_comb; ++i){
        std::vector<int> comb(zero);
        int digit(charge_range.size()-1);
        long n10 = i;
        while (n10 != 0){
            int a = n10 % charge_range[digit].size();
            comb[digit] = a;
            n10 /= charge_range[digit].size();
            --digit;
        }
        std::vector<double> charge_comb;
        int index(0);
        for (int j = 0; j < comb.size(); ++j){
            charge_comb.push_back(candidate[index + comb[j]]);
            index += n_charge_array[j];
        }
        permutation.push_back(charge_comb);
    }

    return permutation;
}

std::vector<double> prepare_charge(const std::vector<double>& charge_candidate,
    const std::vector<int>& spin_array, 
    const std::vector<int>& spin_structure, 
    const std::vector<int>& n_atoms){

    std::vector<double> charge;

    double charge_remaining(0.0);
    for (int i = 0; i < charge_candidate.size(); ++i){
        if (i < spin_array.size()){
            for (int j = 0; j < spin_structure.size(); ++j){
                if (spin_array[i] == spin_structure[j]){
                    charge_remaining -= charge_candidate[i] * n_atoms[j]; 
                    break;
                }
            }

        }
        else {
            charge_remaining -= charge_candidate[i] * n_atoms[i]; 
        }
    }

    std::vector<double> charge_array(charge_candidate);
    charge_array.push_back(charge_remaining / n_atoms[charge_candidate.size()]);

    for (int i = 0; i < spin_structure.size(); ++i){
        for (int j = 0; j < charge_array.size(); ++j){
            if (spin_structure[i] == spin_array[j]){
                for (int k = 0; k < n_atoms[i]; ++k){
                    charge.push_back(charge_array[j]);
                }
                break;
            }
        }
    }
    for (int i = spin_structure.size(); i < charge_array.size(); ++i){
        for (int j = 0; j < n_atoms[i]; ++j){
            charge.push_back(charge_array[i]);
        }
    }

    return charge;
}

int main(){
    
    dvector energy_train, energy_test;
    dmatrix correlation_train, correlation_test, weight_train;
    std::vector<std::string> structure_index_train, structure_index_test;

    // reading data for training
    Input ip_train;
    ip_train.read_correlation("train/correlation");
    ip_train.read_energy("train/dft_energy");
    energy_train = ip_train.get_energy();
    correlation_train = ip_train.get_correlation();
    structure_index_train = ip_train.get_energy_structure_index();
//    weight_train = ublas::identity_matrix<double>(correlation_train.size1());

    // reading data for testing
    Input ip_test;
    ip_test.read_correlation("test/correlation");
    ip_test.read_energy("test/dft_energy");
    energy_test = ip_test.get_energy();
    correlation_test = ip_test.get_correlation();
    structure_index_test = ip_test.get_energy_structure_index();

    // reading opt.in
    const double pi(3.14159265358979323846);

    int n_system, n_atom_prim_first_spin;
    std::vector<int> cluster_index, spin_array;
    vector2_d charge_candidate_array;

    double weight_factor, accuracy, eta, rmax, gmax;

    Parse_input input("opt.in");
    input.assignNeed("cluster_index", cluster_index);
    input.assignNeed("spin", spin_array);
    input.assignNeed("n_charge", n_system);
    input.assignNeed("n_prim_first_spin", n_atom_prim_first_spin);

    for (int i = 0; i < n_system - 1; ++i){
        std::vector<double> charge_input; // init, interval, final
        std::string str1,str2;
        std::stringstream ss;
        str1 = "charge";
        ss << i + 1;
        ss >> str2;
        str1 += str2;
        input.assignNeed(str1.c_str(), charge_input);
        if (charge_input.size() == 3){
            std::vector<double> charge_candidate;
            double tmp(charge_input[0]);
            while (tmp < charge_input[2] + 1e-12){
                charge_candidate.push_back(tmp);
                tmp += charge_input[1];
            }
            charge_candidate_array.push_back(charge_candidate);
        }
        else {
            std::cerr << "  " << str1 << " should have three components."
                << std::endl;
            exit(8);
        }
    }

    input.assign("weight_factor", weight_factor, 0.001);
    input.assign("accuracy", accuracy, 1e-12);
    
    // combinations of effective charges
    vector2_d permutation = get_charge_permutation(charge_candidate_array);

    // ewald calculations and eci estimation
    for (int i = 0; i < permutation.size(); ++i){
        std::vector<double> ewald_energy_array_train, ewald_energy_array_test;
        for (int j = 0; j < structure_index_train.size(); ++j){
            std::string name("train/");
            name += structure_index_train[j];
            name += "/structure_";
            name += structure_index_train[j];

            Parse_structure structure(name.c_str());
            const dmatrix axis = structure.get_axis();
            const std::vector<int> n_atoms = structure.get_num_atoms();
            const std::vector<dvector> position = structure.get_position();
            const double volume = structure.get_volume();
            const dmatrix reciprocal_axis = structure.get_reciprocal_axis();
            double eta = pow ((position.size() * weight_factor
                        * pow(pi,3) / volume), 0.3333333);
            double rmax = sqrt(-log(accuracy) / eta );
            double gmax =  2.0 * sqrt(eta) * sqrt(-log(accuracy));

            std::string name_correlation("train/");
            name_correlation += structure_index_train[j];
            name_correlation += "/correlation.in";
            Parse_input input2(name_correlation.c_str());
            std::vector<int> spin_structure;
            input2.assignNeed("spin", spin_structure);

            std::vector<double> charge = prepare_charge
                (permutation[i], spin_array, spin_structure, n_atoms);

            Ewald ewald
                (axis, position, volume, reciprocal_axis, charge,
                 eta, rmax, gmax);
            ewald.calc_energy();
            double energy = ewald.get_energy();

            double n_cell;
            for (int k = 0; k < spin_structure.size(); ++k){
                if (spin_structure[k] == spin_array[0]){
                    n_cell = double(n_atoms[k]) 
                        / double(n_atom_prim_first_spin);
                }
            }
            ewald_energy_array_train.push_back(energy/n_cell);
            //std::cout.precision(12);
            //std::cout << energy/n_cell << " " << structure_index_train[j] 
            //<< std::endl;
        }

        for (int j = 0; j < structure_index_test.size(); ++j){
            std::string name("test/");
            name += structure_index_test[j];
            name += "/structure_";
            name += structure_index_test[j];

            Parse_structure structure(name.c_str());
            const dmatrix axis = structure.get_axis();
            const std::vector<int> n_atoms = structure.get_num_atoms();
            const std::vector<dvector> position = structure.get_position();
            const double volume = structure.get_volume();
            const dmatrix reciprocal_axis = structure.get_reciprocal_axis();
            double eta = pow ((position.size() * weight_factor
                        * pow(pi,3) / volume), 0.3333333);
            double rmax = sqrt(-log(accuracy) / eta );
            double gmax =  2.0 * sqrt(eta) * sqrt(-log(accuracy));

            std::string name_correlation("test/");
            name_correlation += structure_index_test[j];
            name_correlation += "/correlation.in";
            Parse_input input2(name_correlation.c_str());
            std::vector<int> spin_structure;
            input2.assignNeed("spin", spin_structure);

            std::vector<double> charge = prepare_charge
                (permutation[i], spin_array, spin_structure, n_atoms);

            Ewald ewald
                (axis, position, volume, reciprocal_axis, charge,
                 eta, rmax, gmax);
            ewald.calc_energy();
            double energy = ewald.get_energy();

            double n_cell;
            for (int k = 0; k < spin_structure.size(); ++k){
                if (spin_structure[k] == spin_array[0]){
                    n_cell = double(n_atoms[k]) 
                        / double(n_atom_prim_first_spin);
                }
            }
            ewald_energy_array_test.push_back(energy/n_cell);
        }

        int n_structure_train = correlation_train.size1();
        int n_structure_test = correlation_test.size1();
        int n_cluster = cluster_index.size() + 1;

        dmatrix correlation_calc_train(n_structure_train, n_cluster);
        dmatrix correlation_calc_test(n_structure_test, n_cluster);
        dmatrix weight_calc_train 
            = ublas::identity_matrix<double>(n_structure_train);

        for (int j = 0; j < n_structure_train; ++j){
            for (int k = 0; k < n_cluster-1; ++k){
                correlation_calc_train(j,k) 
                    = correlation_train(j,cluster_index[k]);
            }
            correlation_calc_train(j,n_cluster-1) = ewald_energy_array_train[j];
        }

        for (int j = 0; j < n_structure_test; ++j){
            for (int k = 0; k < n_cluster-1; ++k){
                correlation_calc_test(j,k) 
                    = correlation_test(j,cluster_index[k]);
            }
            correlation_calc_test(j,n_cluster-1) = ewald_energy_array_test[j];
        }

        Least_squares least;
        int check = least.fit(energy_train, correlation_calc_train, 
            weight_calc_train);
        dvector eci;
        double cv_score, error;
        if (check == 1){
            eci = least.get_eci();

            CV cv;
            cv_score = cv.calc_cv_score(energy_train, 
                correlation_calc_train, weight_calc_train, eci);

            dvector pred_energy = ublas::prod(correlation_calc_test, eci);
            dvector diff = pred_energy - energy_test;
            double sq_sum = ublas::inner_prod(diff, diff);
            error = sqrt (sq_sum / double(n_structure_test));

            //std::cout << pred_energy << std::endl;
        }

        else {
            std::cerr << " matrix of correlation functions is singular. "
                << std::endl;
            exit(8);
        }
/*
        double energy_tmp(0.0);
        for (int n = 0; n < 1; ++n){
            for (int j = 0; j < correlation_calc.size2(); ++j){
                energy_tmp += eci(j) * correlation_calc(n, j);
            }
        }
*/
        for (int j = 0; j < permutation[i].size(); ++j){
            std::cout << permutation[i][j] << " ";
        }
        std::cout << std::endl;
        std::cout << eci << std::endl;
        std::cout.precision(14);
        std::cout << " cv_score = " << cv_score << std::endl;
        std::cout << " prediction_error = " << error << std::endl;
//        std::cout << cv_score << std::endl;

    }

    return 0;

}

