/****************************************************************************
  Copyright (c) 2014
  Silke Horn, Andreas Paffenholz (Technische Universitaet Darmstadt, Germany)
  http://www.polymake.org

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2, or (at your option) any
  later version: http://www.gnu.org/licenses/gpl.txt.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
****************************************************************************/



#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Integer.h>

#include <vector>

#include <polymake/common/flint_functions.h>
#include <polymake/polytope/DistanceMatrixPermutation.h>


namespace polymake { 

  namespace polytope {

    // returns a permutation pi for the vector v such that pi(v) 
    // is lex maixmal among all permutations of the entries of v
    std::pair<Vector<Integer>, std::vector<int> > FindMaxPermutation ( const Vector<Integer> & v ) {

      const int dim = v.size();
      Vector<Integer> base(v);
      std::vector<int> perm(dim);
      int k = 0;
      for(std::vector<int>::iterator it = perm.begin(); it != perm.end(); ++it){
        *it = k++;
      }

      for ( int i = 1; i < dim; ++i ) {
	int j = i;
	Integer h = base[i];
	int r = perm[i];
	while ( j > 0 && base[j-1] < h ) {
	  base[j] = base[j-1];
	  perm[j] = perm[j-1];
	  j--;
	}
	base[j] = h;
	perm[j] =r;
      }



/*
      Vector<int> perm(dim);
      for ( int i = 0; i < dim; i++ ) {
	perm[revperm[i]] = i;
      }
*/

      std::pair<Vector<Integer>, std::vector<int> > ret;
      ret.first = base;
      ret.second = perm;

      return ret;
    }


    std::vector<DistanceMatrixPermutation> get_all_permutations_for_row ( int i, const DistanceMatrixPermutation & dmp_in, const Matrix<Integer> & A ) {

      std::vector<DistanceMatrixPermutation> dmp_list;
      Vector<Integer> base(A.cols());

      for ( int j = i; j < A.rows(); ++j ) {
      std::pair<Vector<Integer>, std::vector<int> > perm = FindMaxPermutation(A.row(j));
	if ( perm.first < base ) continue;

	std::vector<int> rperm = dmp_in.get_rperm();
	rperm.push_back(j);
	std::vector<int> blocks;

	DistanceMatrixPermutation dmp(rperm,perm.second,blocks);

	if ( perm.first > base ) {
	  base = perm.first;
	  dmp_list.clear();
	}

	dmp_list.push_back(dmp);
      }

      return dmp_list;
    }



    
    Matrix<Integer> lattice_normal_form ( const perl::Object & p ) {
      
      if ( ! p.give("LATTICE") ) 
	throw std::runtime_error("the given polytope is not a lattice polytope");

      Matrix<Integer> A = p.give("FACET_VERTEX_LATTICE_DISTANCES");

      std::pair<Vector<Integer>,std::vector<int> > fmp = FindMaxPermutation(A.row(0));
      cout << fmp.first << endl;
      std::vector<int> perm = fmp.second;
      for ( int j = 0; j < perm.size(); ++j ) 
	cout << perm[j] << " ";
      cout << endl;
      
      

      std::vector<int> rperm;
      std::vector<int> cperm;
      std::vector<int> blocks;
      DistanceMatrixPermutation dmp(rperm,cperm,blocks);
      
      std::vector<DistanceMatrixPermutation> dmp_list = get_all_permutations_for_row(0,dmp,A);

      for ( int i = 0; i < dmp_list.size(); ++i ) {
	cout << dmp_list[i] << endl;
      }

      return common::HermiteNormalForm(A);
    } 
    
    UserFunction4perl(" ", &lattice_normal_form, "lattice_normal_form(Polytope)");
        
  }
}
