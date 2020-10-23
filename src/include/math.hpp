/****************************************************************************

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

 ******************************************************************************/

#ifndef MATH_HPP
#define MATH_HPP

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

namespace ublas = boost::numeric::ublas;

namespace math {
    template <class M, class MI> 
        int invert(const M& normal, MI& inverse) throw();
    template <class M> double determinant(const M& m);
    template <class V1, class V2>
        ublas::vector<typename ublas::promote_traits<typename V1::value_type, 
        typename V2::value_type>::promote_type>
            cross_prod(const V1& lhs, const V2& rhs);
}

template <class M, class MI>
int math::invert(const M& m, MI& mi) throw()
{
    namespace ublas = boost::numeric::ublas;

    BOOST_UBLAS_CHECK(m.size1() == m.size2(), ublas::external_logic());

    ublas::matrix<double>       lhs(m);
    ublas::matrix<double>       rhs(ublas::identity_matrix<double>(m.size1()));
    ublas::permutation_matrix<> pm(m.size1());

    try {
        ublas::lu_factorize(lhs, pm);
    }
    catch (ublas::internal_logic& e){
        return 0;
    }
    catch (ublas::singular& e){
        return 0;
    }

    try {
        ublas::lu_substitute(lhs, pm, rhs);
    }
    catch (ublas::internal_logic& e){
        return 0;
    }
    catch (ublas::singular& e){
        return 0;
    }

    mi.resize(rhs.size1(), rhs.size2(), false);
    mi.assign(rhs);

    return 1;

}
/*
template <class M, class MI>
int math::invert(const M& normal, MI& inverse) throw()
{

    inverse = normal;

    int signum;
    gsl_permutation *p = gsl_permutation_alloc(normal.size1());
    gsl_matrix *tmp, *dest;

    tmp = gsl_matrix_alloc(normal.size1(), normal.size1());
    dest = gsl_matrix_alloc(normal.size1(), normal.size1());

    for (int i = 0; i < normal.size1(); ++i){
        for (int j = 0; j < normal.size1(); ++j){
            gsl_matrix_set (tmp, i , j, normal(i,j));
        }
    }

    int status;
    status = gsl_linalg_LU_decomp(tmp, p, &signum);
    gsl_linalg_LU_invert(tmp, p, dest);

    gsl_matrix_free(tmp);
    gsl_permutation_free(p);

    for (int i = 0; i < normal.size1(); ++i){
        for (int j = 0; j < normal.size1(); ++j){
            inverse(i,j) = gsl_matrix_get(dest,i,j);
        }
        if (fabs(inverse(i,i)) > 1e10){
            return 0;
        }
    }

    gsl_matrix_free(dest);

    return 1;
}
*/
template <class M>
double math::determinant(const M& m)
{
    namespace ublas = boost::numeric::ublas;

    BOOST_UBLAS_CHECK(m.size1() == m.size2(), ublas::external_logic());

    ublas::matrix<double>       lu(m);
    ublas::permutation_matrix<> pm(m.size1());

    ublas::lu_factorize(lu, pm);

    double det(1);

    typedef ublas::permutation_matrix<>::size_type size_type;

    for (size_type i = 0; i < pm.size(); ++i) {
        det *= (i == pm(i)) ? +lu(i, i) : -lu(i, i);
    }

    return det;
}



template <class V1, class V2>
ublas::vector<typename ublas::promote_traits<typename V1::value_type, 
typename V2::value_type>::promote_type>
math::cross_prod(const V1& lhs, const V2& rhs)
{
    BOOST_UBLAS_CHECK(lhs.size() == 3, ublas::external_logic());
    BOOST_UBLAS_CHECK(rhs.size() == 3, ublas::external_logic());

    typedef typename ublas::promote_traits<typename V1::value_type,
        typename V2::value_type>::promote_type promote_type;
    ublas::vector<promote_type> temporary(3);

    temporary(0) = lhs(1) * rhs(2) - lhs(2) * rhs(1);
    temporary(1) = lhs(2) * rhs(0) - lhs(0) * rhs(2);
    temporary(2) = lhs(0) * rhs(1) - lhs(1) * rhs(0);

    return temporary;
}

#endif

