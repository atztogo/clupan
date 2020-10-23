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


	    Class for calculating cross validation score 

 **************************************************************************/

#include "cv.h"

CV::CV(){}
CV::~CV(){}

double CV::calc_cv_score(const dvector& energy, 
        const dmatrix& correlation, const dmatrix& weight, 
        const dvector& eci){

    dmatrix prod_mat(correlation);
    for (int i = 0; i < correlation.size1(); ++i){
        for (int j = 0; j < correlation.size2(); ++j){
            prod_mat(i,j) *= weight(i,i);
        }
    }

    dmatrix trans_mat = ublas::trans(correlation);
    dmatrix prod_mat2 = ublas::prod(trans_mat, prod_mat);
    dmatrix inverse;
    int check = math::invert(prod_mat2, inverse);

    int n_struct = correlation.size1();
    int n_cluster = correlation.size2();
    double denom;
    if (check == 1){
        double square_cv = 0.0;
        for (int i = 0; i < n_struct; ++i){
            double sqrt_weight = sqrt(weight(i,i));
            dvector correlation_i = ublas::row(correlation,i);
            dvector tmp_vec = sqrt_weight 
                * ublas::prod(correlation_i, inverse);
            dvector tmp_vec2 = ublas::trans(correlation_i) * sqrt_weight;
            long double denom = 1 - inner_prod(tmp_vec,tmp_vec2);
            if (denom < 1e-13 or fabs(denom) > 1){
                goto end;
            }
            long double predict_energy(0.0);
            for (int j = 0; j < n_cluster; ++j){
                predict_energy += eci(j) * correlation(i,j);
            }
            long double tmp = (energy(i) - predict_energy) / denom;
            square_cv += weight(i,i) * pow(tmp, 2);
        }
        cv_score = sqrt(square_cv / n_struct);
    }
    else {
        end:;
        cv_score = 1e10;
    }


    return cv_score;

}

std::vector<double> CV::calc_cv_casp(const dvector& energy, 
        const dmatrix& correlation, const dvector& eci, 
        const std::vector<int>& group){

    int max_group = *std::max_element(group.begin(), group.end());
    for (int i = 0; i < max_group; ++i){
        cv_casp.push_back(0.0);
        n_structure_in_group.push_back(0);
    }

    dmatrix trans_mat = ublas::trans(correlation);
    dmatrix prod_mat2 = ublas::prod(trans_mat, correlation);
    dmatrix inverse;
    int check = math::invert(prod_mat2, inverse);

    int n_struct = correlation.size1();
    int n_cluster = correlation.size2();
    if (check == 1){
        double square_cv = 0.0;
        for (int i = 0; i < n_struct; ++i){
            dvector correlation_i = ublas::row(correlation,i);
            dvector tmp_vec = ublas::prod(correlation_i, inverse);
            dvector tmp_vec2 = ublas::trans(correlation_i);
            double denom = 1 - inner_prod(tmp_vec,tmp_vec2);
            double predict_energy(0.0);
            for (int j = 0; j < n_cluster; ++j){
                predict_energy += eci(j) * correlation(i,j);
            }
            double tmp = (energy(i) - predict_energy) / denom;
            cv_casp[group[i]-1] += pow(tmp, 2);
            n_structure_in_group[group[i]-1] += 1;
        }
        for (int i = 0; i < cv_casp.size(); ++i){
            cv_casp[i] = sqrt(cv_casp[i] / n_structure_in_group[i]);
        }
    }
    else {
        std::cerr << " CV-CASP cannot be obtained. " << std::endl;
        exit(8);
    }

    return cv_casp;

}

const std::vector<int>& CV::get_n_structure_in_group() const{
    return n_structure_in_group;
}
/*
std::vector<double> CV::calc_cv_casp2(const dvector& energy, 
        const dmatrix& correlation, const std::vector<int>& group){

    int max_group = *std::max_element(group.begin(), group.end());
    std::vector<double> cv_casp;
    std::vector<int> n_structure_in_group;
    for (int i = 0; i < max_group; ++i){
        cv_casp.push_back(0.0);
        n_structure_in_group.push_back(0);
    }

    int n_struct = correlation.size1();
    int n_cluster = correlation.size2();

    for (int n = 0; n < n_struct; ++n){
        dmatrix correlation_tmp(n_struct - 1, n_cluster);
        dvector energy_tmp(n_struct-1);
        dmatrix weight_tmp = ublas::identity_matrix<double> (n_struct-1);
        int ncount(0);
        for (int i = 0; i < n_struct; ++i){
            if (n != i){
                for (int j = 0; j < n_cluster; ++j){
                    correlation_tmp(ncount,j) = correlation(i,j);
                }
                energy_tmp(ncount) = energy(i);
                ++ncount;
            }
        }

        Least_squares least;
        int check = least.fit(energy_tmp, correlation_tmp, weight_tmp);
        dvector eci = least.get_eci();
        dvector correlation_n = ublas::row(correlation,n);
        double energy_ce = ublas::inner_prod(correlation_n, eci);
        double diff_energy = energy_ce - energy(n);
        cv_casp[group[n]-1] += pow(diff_energy, 2);
        n_structure_in_group[group[n]-1] += 1;
    }

    for (int i = 0; i < cv_casp.size(); ++i){
        cv_casp[i] = sqrt(cv_casp[i] / n_structure_in_group[i]);
    }

    return cv_casp;

}
*/
