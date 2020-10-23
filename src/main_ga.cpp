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

	    Main program for genetic algorithm
	
*****************************************************************************/

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_rng.h>

#include "input.h"
#include "output.h"
#include "ga.h"
#include "combination.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> dvector;
typedef ublas::matrix<int> imatrix;
typedef ublas::matrix<double> dmatrix;

typedef std::vector<std::vector<int> > vector2_i;

const double k_B = 0.86171 * pow (10,-4);

int main(){

    Input ip;
    ip.input_ga();

    dvector energy = ip.get_energy();
    dmatrix correlation_all = ip.get_correlation();
    dmatrix weight = ip.get_weight();

    // input parameters
    int nloop = ip.get_nloop();
	std::vector<double> temp_array = ip.get_temp_array();
    int n_cluster = ip.get_n_cluster();
	std::vector<int> base_array = ip.get_base_array();
	std::vector<int> cluster_candidate = ip.get_cluster_candidate();

    int n_pop = ip.get_n_pop();
    int max_iteration = ip.get_max_iteration();
    int n_elite = ip.get_n_elite();
    double p_mating = ip.get_p_mating();
    double p_mutation = ip.get_p_mutation();

    // initialization for generating random numbers
    unsigned long init = time(NULL) + getpid();
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, init);
    
    // random preparation of population 
    std::cout << "  generating population " << std::endl;
    GA ga;
    vector2_i pop_array(n_pop);
    std::vector<double> score_array(n_pop);

	for (int i = 0; i < n_pop; ++i){
        std::vector<int> pop_init;
        double score_init(1e100);
        while (score_init > 1e99 ){
            pop_init = ga.initial_combination(cluster_candidate, n_cluster, r);
            score_init = ga.calc_score(energy, correlation_all, weight, 
                base_array, pop_init);
        }
        pop_array[i] = pop_init;
        score_array[i] = score_init;
    }

    // optimization of population using simulated annealing
    std::cout << "  generating population using simulated annealing" 
        << std::endl;
    if (nloop > 0){
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int i = 0; i < pop_array.size(); ++i){
            std::vector<int> pop_old(pop_array[i]), pop_new;
            double score_old(score_array[i]), score_new, diff_score, prob;
            for (int itemp = 0; itemp < temp_array.size(); ++itemp){
                for (int iloop = 0; iloop < nloop; ++iloop){
                    pop_new = ga.change_combination
                        (cluster_candidate, pop_old, r);
                    score_new = ga.calc_score(energy, correlation_all, 
                            weight, base_array, pop_new);
                    diff_score = score_new - score_old;
                    prob = exp (- diff_score / k_B / temp_array[itemp]);
                    if (prob > gsl_rng_uniform(r)){
                        pop_old = pop_new;
                        score_old = score_new;
                    }
                }
            }
            pop_array[i] = pop_old;
            score_array[i] = score_old;
        }
    }

    std::cout << "  iteration 0 : min_cv_score = " 
        << *std::min_element(score_array.begin(), score_array.end()) 
        << std::endl;

    // output of initial population
    Output op;
    op.output_ga(pop_array, score_array, -1);

    // combination of parents for making child
    Combination comb;
    imatrix comb_pop = comb.get_combination(n_pop, 2);
    int n_pair = comb_pop.size1();

    vector2_i pop_array_old(pop_array);
    std::vector<double> score_array_old(score_array);
    for  (int i_iter = 0; i_iter < max_iteration; ++i_iter){
        vector2_i pop_array_new;
        std::vector<double> score_array_new;

        // making child and calculating cv score for the child
        while (pop_array_new.size() != n_pop){
            int comb1 = gsl_rng_uniform_int(r, n_pair);
            int parent1(comb_pop(comb1,0)), parent2(comb_pop(comb1,1));
            std::vector<int> child = ga.generate_child(pop_array_old[parent1], 
                    pop_array_old[parent2], cluster_candidate, p_mating, 
                    p_mutation, n_elite, r);
            double score_child = ga.calc_score(energy, correlation_all, 
                    weight, base_array, child);
            if (score_child <= score_array_old[parent1] 
                    and score_child <= score_array_old[parent2]){
                pop_array_new.push_back(child);
                score_array_new.push_back(score_child);
            }
        }

        // combining and sorting old population and new population
        std::multimap<double, std::vector<int> > pop_array_sort;
        std::pair<double, std::vector<int> > tmp_pair;
        for (int i = 0; i < pop_array_new.size(); ++i){
            tmp_pair.first = score_array_new[i];
            tmp_pair.second = pop_array_new[i];
            pop_array_sort.insert(tmp_pair);
        }
        for (int i = 0; i < pop_array_old.size(); ++i){
            tmp_pair.first = score_array_old[i];
            tmp_pair.second = pop_array_old[i];
            pop_array_sort.insert(tmp_pair);
        }

        std::multimap<double, std::vector<int> >::iterator it_map 
            = pop_array_sort.begin();

        // elite selection
        vector2_i pop_array_elite;
        std::vector<double> score_array_elite;
        for (int i = 0; i < n_elite; ++i){
            score_array_elite.push_back((*it_map).first);
            pop_array_elite.push_back((*it_map).second);
            pop_array_sort.erase(it_map);
            ++it_map;
        }

        for (int i = 0; i < n_pop - n_elite; ++i){
            int rand1 = gsl_rng_uniform_int(r, pop_array_sort.size());
            it_map = pop_array_sort.begin();
            for (int j = 0; j < rand1; ++j){
                ++it_map;
            }
            score_array_elite.push_back((*it_map).first);
            pop_array_elite.push_back((*it_map).second);
            pop_array_sort.erase(it_map);
        }

        pop_array_old = pop_array_elite;
        score_array_old = score_array_elite;

        std::cout << "  iteration " << i_iter + 1 
            << " : min_cv_score = " << score_array_old[0] << std::endl;

        Output op;
        op.output_ga(pop_array_old, score_array_old, i_iter);
    }

    return 0;

}

