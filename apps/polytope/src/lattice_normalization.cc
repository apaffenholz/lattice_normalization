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

#define DEBUG

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Integer.h>

#include <vector>

#include <polymake/common/flint_functions.h>
#include <polymake/polytope/DistanceMatrixPermutation.h>


namespace polymake { 

  namespace polytope {

    // apply a permutation perm to a vector v
    Vector<Integer> apply_permutation ( const Vector<Integer> & v, std::vector<int> perm ) {
      
      Vector<Integer> w(v.size());
      for ( int i = 0; i < v.size(); ++i )
	w[i] = v[perm[i]];
      
      return w;
    }


    // returns a permutation pi for the vector v such that pi(v) 
    // is lex maixmal among all permutations of the entries of v
    std::pair<Vector<Integer>, std::vector<int> > FindMaxPermutation ( const Vector<Integer> & v,
								       const std::vector<int> perm_in,
								       const std::vector<int> blocks ) {


#ifdef DEBUG
      cout << "[FindMaxPermutation] incoming data: " << v << endl << "permutation: " << endl;
      for ( int i = 0; i < perm_in.size(); i++ ) {
	cout << perm_in[i] << " ";
      }
      cout << endl << "blocks: " << endl;
      for ( int i = 0; i < blocks.size(); i++ ) {
	cout << blocks[i] << " ";
      }
      cout << endl << "-----------------------" << endl;
#endif

      
      const int dim = v.size();
      Vector<Integer> base = apply_permutation(v, perm_in);
      std::vector<int> tperm(dim);
      int k = 0;
      for(std::vector<int>::iterator it = tperm.begin(); it != tperm.end(); ++it){
        *it = k++;
      }

      int j = 0;
      for ( int i = 0; i < blocks.size(); ++i ) {
	int k = j;
	j++;
	for ( ; j <= blocks[i]; ++j ) {
	  int l = j;
	  Integer h = base[j];
	  int r = tperm[j];
	  while ( l > k && base[l-1] < h ) {
	    base[l] = base[l-1];
	    tperm[l] = tperm[l-1];
	    l--;
	  }
	  base[l] = h;
	  tperm[l] =r;
	}
      }


      std::vector<int> perm(dim);
      for ( int i = 0; i < dim; i++ ) {
	perm[i] = perm_in[tperm[i]];
      }

      std::pair<Vector<Integer>, std::vector<int> > ret;
      ret.first = base;
      ret.second = perm;

#ifdef DEBUG
      cout << "[FindMaxPermutation] returning base: " << base << endl << "permutation: " << endl;
      for ( int i = 0; i < dim; i++ ) {
	cout << perm[i] << " ";
      }
      cout << endl << "-----------------------" << endl;
#endif


      return ret;
    }


    // given an initial permutation dmp_in with partially filled rperm
    // compute all possible new permutations for rows with index i or higher
    // FIXME: needs to be rewritten as the function does not respect the initial row permutation in dmp_in
    std::vector<DistanceMatrixPermutation> get_all_permutations_for_row ( int i, const DistanceMatrixPermutation & dmp_in, const Matrix<Integer> & A ) {

      std::vector<DistanceMatrixPermutation> dmp_list;
      Vector<Integer> base(A.cols());

      for ( int j = i; j < A.rows(); ++j ) {
	std::pair<Vector<Integer>, std::vector<int> > perm = FindMaxPermutation(A.row(j), dmp_in.get_cperm(), dmp_in.get_blocks() );
	if ( perm.first < base ) continue;

	std::vector<int> rperm = dmp_in.get_rperm();
	rperm.push_back(j);
	std::vector<int> blocks;

	Integer e = perm.first[0];
	for( int k = 1; k < A.cols(); ++k )
	  if ( perm.first[k] == e ) {
	    continue;
	  } else {
	    blocks.push_back(k-1);
	    e = perm.first[k];
	  }
	blocks.push_back(A.cols()-1);

	DistanceMatrixPermutation dmp(rperm,perm.second,blocks);

	if ( perm.first > base ) {
	  base = perm.first;
	  dmp_list.clear();
	}

	dmp_list.push_back(dmp);
      }

      return dmp_list;
    }

    
    // return the lattice normal form of a lattice polytope
    Matrix<Integer> lattice_normal_form ( const perl::Object & p ) {

      // polytope must be lattice, otherwise the distance matrix is not defined
      if ( ! p.give("LATTICE") ) 
	throw std::runtime_error("the given polytope is not a lattice polytope");

      Matrix<Integer> A = p.give("FACET_VERTEX_LATTICE_DISTANCES");

      // initialize a matrix permutation with the identity permutation 
      std::vector<int> rperm;
      std::vector<int> cperm(A.cols());
      int k = 0;
      for(std::vector<int>::iterator it = cperm.begin(); it != cperm.end(); ++it){
        *it = k++;
      }

      std::vector<int> blocks(1);
      blocks[0] = A.cols()-1;
      
      DistanceMatrixPermutation dmp(rperm,cperm,blocks);

      
#ifdef DEBUG
      cout << "[lattice_normalization] initial permutation matrix: " << endl << dmp << endl << "-----------------------" << endl;
#endif

      // try for row 0
      std::vector<DistanceMatrixPermutation> dmp_list = get_all_permutations_for_row(0,dmp,A);

      cout << "current base vector for row 0: " << apply_permutation(A.row(0),dmp_list[0].get_cperm()) << endl;
      for ( int i = 0; i < dmp_list.size(); ++i ) {
	cout << dmp_list[i] << endl;
      }

      // try for row 1
      // currently other rows won't work as we don't yet respect rperm in choosing the next row
      std::vector<DistanceMatrixPermutation> dmp_list_1 = get_all_permutations_for_row(1,dmp_list[0],A);

      cout << "current base vector for row 1: " << apply_permutation(A.row(1),dmp_list_1[0].get_cperm()) << endl;
      for ( int i = 0; i < dmp_list_1.size(); ++i ) {
	cout << dmp_list_1[i] << endl;
      }
      
      // just to return something
      return common::flint::HermiteNormalForm(A);
    } 
    
    UserFunction4perl(" ", &lattice_normal_form, "lattice_normal_form(Polytope)");
        
  }
}
