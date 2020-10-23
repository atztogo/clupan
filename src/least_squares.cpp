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


	    Class for least squares fitting

 ****************************************************************************/

#include "least_squares.h"

Least_squares::Least_squares(){}

Least_squares::~Least_squares(){}
/*
int Least_squares::fit(const dvector& energy, 
    const dmatrix& correlation, const dmatrix& weight){

    dmatrix trans_correlation = ublas::trans(correlation);
    dmatrix prod_mat = ublas::prod(trans_correlation, weight);
    dmatrix prod_mat2 = ublas::prod(prod_mat, correlation);

    dmatrix inverse_prod_mat2;
    int check = math::invert(prod_mat2, inverse_prod_mat2);
    if (check == 1){
        dmatrix prod_mat3 = ublas::prod(inverse_prod_mat2, prod_mat);
        eci = ublas::prod(prod_mat3, energy);
        dvector pred_energy = ublas::prod(correlation, eci);
        dvector diff = energy - pred_energy;
        dvector weight_diff(diff);
        for (int i = 0; i < diff.size(); ++i){
            weight_diff(i) = weight(i,i) * diff(i);
        }
        square_error = ublas::inner_prod(diff, weight_diff);
    }

    return check;
}
*/

int Least_squares::fit(const dvector& energy, 
        const dmatrix& correlation, const dmatrix& weight){

    int check(1);

    const int n_struct(correlation.size1()), n_cluster(correlation.size2());
    double chisq;

    //dmatrix covtemp(nCluster, nCluster);

    gsl_matrix *correlation_gsl, *cov_gsl;
    gsl_vector *energy_gsl, *weight_gsl, *eci_gsl;

    correlation_gsl = gsl_matrix_alloc (n_struct, n_cluster);
    energy_gsl = gsl_vector_alloc (n_struct);
    weight_gsl = gsl_vector_alloc (n_struct);

    eci_gsl = gsl_vector_alloc (n_cluster);
    cov_gsl = gsl_matrix_alloc (n_cluster, n_cluster);

    for (int i = 0; i < n_struct; ++i){
        gsl_vector_set (energy_gsl, i, energy(i));
        gsl_vector_set (weight_gsl, i, weight(i,i));
        for (int j = 0; j < n_cluster; ++j){
            gsl_matrix_set (correlation_gsl, i , j,
                    correlation(i,j));
        }
    }

    gsl_multifit_linear_workspace * work
        = gsl_multifit_linear_alloc (n_struct, n_cluster);

    int status;
    gsl_set_error_handler_off();

    status = gsl_multifit_wlinear (correlation_gsl, weight_gsl,
            energy_gsl, eci_gsl, cov_gsl, &chisq, work);
    
    if (status == GSL_EINVAL) {
        std::cerr << " Warning : matrix of correlation functions is singular. "
            << std::endl;
        //fprintf(stderr, "failed, gsl_errno=%d\n", status);
    } 
    else {
        //fprintf(stderr, "failed, gsl_errno=%d\n", status);
    }
    //exit (-1);
 
    gsl_multifit_linear_free (work);

    eci.resize(n_cluster);
    for (int i = 0; i < n_cluster; ++i)
        eci(i) = gsl_vector_get(eci_gsl,(i));

    square_error = chisq;

    gsl_matrix_free (correlation_gsl);
    gsl_vector_free (energy_gsl);
    gsl_vector_free (weight_gsl);
    gsl_matrix_free (cov_gsl);
    gsl_vector_free (eci_gsl);

    return check;
}

const dvector& Least_squares::get_eci() const{

    return eci;
}

const double& Least_squares::get_square_error(){
    
    return square_error;
}

