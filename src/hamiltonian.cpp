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

	    Class for constructing hamiltonian
	
****************************************************************************/

#include "hamiltonian.h"

Hamiltonian::Hamiltonian(const std::vector<int>& spin){

	onsiteFunction = gramSchmidt(spin);

}

Hamiltonian::~Hamiltonian(){}

dmatrix Hamiltonian::gramSchmidt(const std::vector<int>& spin){

	// cluster approach to order-disorder transformations p.50

	// startFunction : coefficient of w_k : 1, sigma, sigma^2 ...
	// basis : coefficient of g_k
	// localFunction : coefficient of phi_k

	int nsize = spin.size();
	const dmatrix startFunction = ublas::identity_matrix<double> (nsize);

	dmatrix localFunction(ublas::zero_matrix<double> (nsize,nsize));

	for (int i = 0; i < nsize; ++i){
		dvector basis = ublas::row(startFunction,i);
		for (int j = 0; j < i; ++j){
			basis -= innerProduct(spin, 
					ublas::row(localFunction,j), 
					ublas::row(startFunction,i)) 
				* ublas::row(localFunction,j);
		}
		double norm = sqrt(innerProduct(spin, basis, basis));
		for (int j = 0; j < nsize; ++j){
			localFunction(i, j) = basis(j) / norm;
		}
	}

	return localFunction;	
}

// local function of calculating inner product
dvector Hamiltonian::functionProduct(const dvector& func1, 
					const dvector& func2){

	dvector funcProd
		(ublas::zero_vector<double>(func1.size() + func2.size() - 1));
	
	for (int i = 0; i < func1.size(); ++i){
		for (int j = 0; j < func2.size(); ++j)
			funcProd(i+j) += func1(i) * func2(j);
	}

	return funcProd;
}

// local function of calculating inner product
double Hamiltonian::innerProduct(const std::vector<int>& spin, 
		const dvector& func1, const dvector& func2){

	const int normal = spin.size();
	const dvector funcProd(functionProduct(func1, func2));

	double inner(0);
	for (int i = 0; i < normal; ++i){
		for (int j = 0; j < funcProd.size(); ++j)
			inner += funcProd(j) * pow(spin[i],j); 
	}

	inner /= double(normal);

	return inner;
	
}

const dmatrix& Hamiltonian::get_onsite_function() const{
	
	return onsiteFunction;
}
