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

	
	    Class for genetic algorithm

******************************************************************************/

#include "ga.h"

GA::GA(){}

GA::~GA(){}

std::vector<int> GA::initial_combination
(const std::vector<int>& candidate, const int& n_cluster, gsl_rng * r){

    std::vector<int> combination(n_cluster);

    std::set<int> combination_set;
    while (combination_set.size() < n_cluster){
        int rand1 = gsl_rng_uniform_int(r, candidate.size());
        combination_set.insert(candidate[rand1]);
    }

    std::copy(combination_set.begin(), combination_set.end(), 
        combination.begin());

    return combination;
}

double GA::calc_score(const dvector& energy, 
        const dmatrix& correlation_all, 
        const dmatrix& weight, const std::vector<int>& base_array, 
        const std::vector<int>& combination){

    double cv_score;

    const int n_struct(correlation_all.size1());
    dmatrix correlation(n_struct, base_array.size() + combination.size());
    for (int i = 0; i < n_struct; ++i){
        for (int j = 0; j < base_array.size(); ++j){
            correlation(i,j) = correlation_all(i, base_array[j]);
        }
        for (int j = 0; j < combination.size(); ++j){
            correlation(i, base_array.size() + j)
                = correlation_all(i,combination[j]);
        }
    }

    Least_squares least;
    int check = least.fit(energy, correlation, weight);
    if (check == 1){
        dvector eci = least.get_eci();
        CV cv;
        cv_score = cv.calc_cv_score(energy, correlation, weight, eci);
    }
    else {
        cv_score = 1e10;
    }

    return cv_score;
}

std::vector<int> GA::change_combination(
    const std::vector<int>& candidate, 
    const std::vector<int>& combination_old, gsl_rng * r){

    std::vector<int> combination_new;

    int rand1 = gsl_rng_uniform_int(r, combination_old.size());
    combination_new = combination_old;
    combination_new.erase(combination_new.begin() + rand1);
    int rand2 = gsl_rng_uniform_int(r, candidate.size());
    while (find(combination_new.begin(), combination_new.end(), 
        candidate[rand2]) != combination_new.end()) {
        rand2 = gsl_rng_uniform_int(r, candidate.size());
    }
    combination_new.push_back(candidate[rand2]);

    return combination_new;

}

std::vector<int> GA::generate_child
(const std::vector<int>& parent1, const std::vector<int>& parent2,
 const std::vector<int>& candidate, const double& p_mating, 
 const double& p_mutation, const int& n_elite, gsl_rng * r){

    std::vector<int> child;
    const int n_select(parent1.size());

    /// mating
    for  (int i = 0; i < parent1.size(); ++i){
        if (find(parent2.begin(), parent2.end(), parent1[i]) != parent2.end())
            child.push_back(parent1[i]);
        else if (p_mating > gsl_rng_uniform(r))
            child.push_back(parent1[i]);
    }
    for  (int i = 0; i < parent2.size(); ++i){
        if (find(parent1.begin(), parent1.end(), parent2[i]) == parent1.end() 
            and p_mating < gsl_rng_uniform(r)){
            child.push_back(parent2[i]);
        }
    }

    /// mutation
    for  (int i = 0; i < candidate.size(); ++i){
        if (p_mutation > gsl_rng_uniform(r)){
            std::vector<int>::iterator it_vec = 
                find(child.begin(), child.end(), candidate[i]);
            if (it_vec != child.end())
                child.erase(it_vec);
            else 
                child.push_back(candidate[i]);
        }
    }

    /// reduction or addition
    int diff_size = child.size() - n_select;
    if (diff_size > 0){
        while (child.size() != n_select){
            int rand1 = gsl_rng_uniform_int(r, child.size());
            child.erase(child.begin() + rand1);
        }
    }
    else if (diff_size < 0){
        while (child.size() != n_select){
            int rand1 = gsl_rng_uniform_int(r, candidate.size());
            if (find(child.begin(), child.end(), candidate[rand1]) 
                == child.end()){
                child.push_back(candidate[rand1]);
            }
        }
    }

    return child;

}
