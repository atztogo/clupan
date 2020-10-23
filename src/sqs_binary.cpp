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
        
        Class for finding sqs in binary systems

***************************************************************************/

#include "sqs_binary.h"

SQS_binary::SQS_binary(const std::vector<int>& spin_input,
    const dvector& correlation_array_input, 
    const dvector& n_cluster_array_input, 
    std::vector<imatrix>& permutation_array,
    const int& n_step_ann_input, 
    const std::vector<double>& temp_array_input,
    const double& disorder_cf_input, 
    const std::string& criterion_input):
    spin(spin_input), correlation_array(correlation_array_input),
    n_cluster_array(n_cluster_array_input), n_step_ann(n_step_ann_input), 
    temp_array(temp_array_input), disorder_cf(disorder_cf_input),
    criterion(criterion_input){

    // coordination
    set_atomic_coordination(permutation_array);

}

SQS_binary::~SQS_binary(){}

void SQS_binary::cmc(){

    const double k_b(0.00008617343);
    unsigned long init = time(NULL) + getpid();
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, init);

    int n_cluster = correlation_array.size();
    int n_atom_simulation = spin.size();

    std::vector<int> spin_min;
    dvector correlation_array_min;
    double score(0.0), score_min(1e100);
    if (criterion == "rms"){
        for (int i = 0; i < n_cluster; ++i){
            score += pow((correlation_array(i) - disorder_cf), 2.0);
        }
    }
    else if (criterion == "abs"){
        for (int i = 0; i < n_cluster; ++i){
            score += fabs(correlation_array(i) - disorder_cf);
        }
    }
 
    for (int i = 0; i < temp_array.size(); ++i){
        double temperature;
        int n_step(n_step_ann);
        temperature = temp_array[i];
        std::cout << "  temperature = " << temperature << std::endl;
        std::cout << "  n_step = " << n_step << std::endl;
        
        int step_sum(0), atom1, atom2, spin1, spin2;
        double tmp_prob, prob, cut, score_diff;
        dvector diff_spin, spin_old(n_cluster), spin_new(n_cluster), 
            diff_correlation(n_cluster);
        std::vector<int>::iterator it1, it2;
        while (step_sum < n_step){
            atom1 = gsl_rng_uniform_int(r, n_atom_simulation);
            atom2 = gsl_rng_uniform_int(r, n_atom_simulation);
            spin1 = spin[atom1];
            spin2 = spin[atom2];
            if (spin1 != spin2){
                for (int j = 0; j < n_cluster; ++j){
                    spin_old(j) = calc_spin_binary
                        (spin, coordination[atom1][j], 
                        n_cluster_coordination[atom1][j])
                        + calc_spin_binary
                        (spin, coordination[atom2][j], 
                        n_cluster_coordination[atom2][j]);
                }
                
                it1 = spin.begin() + atom1;
                it2 = spin.begin() + atom2;
                iter_swap(it1, it2);
                for (int j = 0; j < n_cluster; ++j){
                    spin_new(j) = calc_spin_binary
                        (spin, coordination[atom1][j], 
                        n_cluster_coordination[atom1][j])
                        + calc_spin_binary
                        (spin, coordination[atom2][j], 
                        n_cluster_coordination[atom2][j]);
                }
                diff_spin = spin_new - spin_old;
                dvector diff_correlation(diff_spin);
                for (int j = 0; j < n_cluster; ++j){
                    diff_correlation(j) /= n_cluster_array(j);
                }

                dvector correlation_array_new 
                    = correlation_array + diff_correlation;
                double score_new(0.0);
                if (criterion == "rms"){
                    for (int j = 0; j < n_cluster; ++j){
                        score_new += pow((correlation_array_new(j) 
                                    - disorder_cf), 2.0);
                    }
                }
                else if (criterion == "abs"){
                    for (int j = 0; j < n_cluster; ++j){
                        score_new 
                            += fabs(correlation_array_new(j) - disorder_cf);
                    }
                }
                score_diff = score_new - score;
                tmp_prob = - score_diff / (k_b * temperature);
                prob = exp(tmp_prob);
                cut =  gsl_rng_uniform(r);
                if (prob > cut){
                    correlation_array += diff_correlation;
                    score = score_new;
                    if (score < score_min){
                        score_min = score;
                        spin_min = spin;
                        correlation_array_min = correlation_array;
                    }
                }
                else {
                    it1 = spin.begin() + atom1;
                    it2 = spin.begin() + atom2;
                    iter_swap(it1, it2);
                }
                ++step_sum;
            }
        }
        if (criterion == "rms"){
            std::cout << "   rms (sqs) = " 
                << sqrt(score/double(n_cluster)) << std::endl;
        }
        else if (criterion == "abs"){
            std::cout << "   mean of abs (sqs) = "
                << score/double(n_cluster) << std::endl;
        }
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

void SQS_binary::set_atomic_coordination(
    std::vector<imatrix>& permutation_array){

    // coordination
    int n_atom_simulation = spin.size();
    coordination.resize(n_atom_simulation);
    n_cluster_coordination.resize(n_atom_simulation);
    for (int i = 0; i < n_atom_simulation; ++i){
        coordination[i].resize(permutation_array.size());
        n_cluster_coordination[i].resize(permutation_array.size());
    }

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int i = 0; i < permutation_array.size(); ++i){
        std::multiset<std::vector<int> > cluster_array;
        for (int j = 0; j < permutation_array[i].size1(); ++j){
            std::multiset<int> cluster_set;
            for (int k = 0; k < permutation_array[i].size2(); ++k){
                cluster_set.insert(permutation_array[i](j,k));
            }
            std::vector<int> cluster(cluster_set.size());
            std::copy(cluster_set.begin(), cluster_set.end(), cluster.begin());
            cluster_array.insert(cluster);
        }
        std::multiset<std::vector<int> >::iterator it = cluster_array.begin();
        while (it != cluster_array.end()){
            std::vector<int> cluster = *it;
            int n_same = cluster_array.count(*it);
            std::set<int> cluster_duplication_check;
            for (int j = 0; j < cluster.size(); ++j){
                cluster_duplication_check.insert(cluster[j]);
            }
            std::set<int>::iterator it_set = cluster_duplication_check.begin();
            while (it_set != cluster_duplication_check.end()){
                coordination[*it_set][i].push_back(cluster);
                n_cluster_coordination[*it_set][i].push_back(n_same);
                ++it_set;
            }
            for (int j = 0; j < n_same; ++j){
                ++it;
            }
        }
    }
}

double SQS_binary::calc_spin_binary(const std::vector<int>& spin,
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

const std::vector<int>& SQS_binary::get_spin() const{
    return spin;
}
const dvector& SQS_binary::get_correlation_array() const{
    return correlation_array;
}
