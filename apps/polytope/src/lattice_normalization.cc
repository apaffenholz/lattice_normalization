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

//#define DEBUG

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Integer.h>

#include <vector>

#include <polymake/common/flint_functions.h>
#include <polymake/polytope/DistanceMatrixPermutation.h>


namespace polymake { 

  namespace polytope {

    // apply a permutation perm to a vector v
    // note that permutations are noted in inverted form: perm[i] is the element that should be placed in position i
    Vector<Integer> apply_permutation ( const Vector<Integer> & v, std::vector<int> perm ) {
      
      Vector<Integer> w(v.size());
      for ( int i = 0; i < v.size(); ++i )
	w[i] = v[perm[i]];
      
      return w;
    }

    

    // returns a permutation pi for the vector v such that pi(v) 
    // is lex maximal among all permutations of the entries of v
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
#endif

      const int dim = v.size();

      // apply the current column permutation to v
      // we can then permute entries of v within each block to generate a lex larger vector
      Vector<Integer> base = apply_permutation(v, perm_in);

      // initial array to store the permutation applied to the permuted row v
      // initialize with identity
      std::vector<int> tperm(dim);
      int k = 0;
      for(std::vector<int>::iterator it = tperm.begin(); it != tperm.end(); ++it)
        *it = k++;

      //find a perm maximising the vector using simple insertion sort
      // FIXME maybe we can do this more efficiently
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

      // concatenate initial permutation with the one computed for the permuted row
      std::vector<int> perm(dim);
      for ( int i = 0; i < dim; i++ ) {
	perm[i] = perm_in[tperm[i]];
      }

      // prepare return
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
    void get_all_permutations_for_row ( int i,
					const DistanceMatrixPermutation & dmp_in,
					const Matrix<Integer> & A,
					std::vector<DistanceMatrixPermutation>& dmp_list,
					Vector<Integer>& base
					) {

      for ( Entire<Set<int> >::const_iterator sit = entire(dmp_in.get_still_available_rows()); !sit.at_end(); ++sit ) {
	std::pair<Vector<Integer>, std::vector<int> > perm = FindMaxPermutation(A.row(*sit), dmp_in.get_cperm(), dmp_in.get_blocks() );

#ifdef DEBUG
	cout << "[get_all_permutations_for_row] found permuted vector: " << perm.first << endl << "comparing to base: " << base << endl;
#endif

	// compare the max permutation achieved from permuting row *sit to the overall max seen so far
	// if to small do nothing
	if ( perm.first < base ) continue;

	// if larger then dicard the previous permutations
	if ( perm.first > base ) {
	  base = perm.first;
	  dmp_list.clear();
	}

	// compute the new block decomposition of the columns
	std::vector<int> rperm = dmp_in.get_rperm();
	rperm.push_back(*sit);
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

	// initialize a new DistancematrixPermutation and add it
	DistanceMatrixPermutation dmp(rperm,perm.second,blocks,dmp_in.get_still_available_rows()-*sit);
	dmp_list.push_back(dmp);
      }

    }

    
    // return the lattice normal form of a lattice polytope
    Matrix<Integer> affine_lattice_normal_form ( const perl::Object & p ) {

      // polytope must be lattice, otherwise the distance matrix is not defined
      if ( ! p.give("LATTICE") ) 
	throw std::runtime_error("the given polytope is not a lattice polytope");

      // get the lattice distances between vertices and facets
      Matrix<Integer> A = p.give("FACET_VERTEX_LATTICE_DISTANCES");

      // get the vertices and dehomogenize
      Matrix<Integer> V = p.give("VERTICES");
      V = V.minor(All,~scalar2set(0));

      // initialize a permutation fo the rows, will be built up sequentially
      std::vector<int> rperm;
      
      // initialize a permutation of the columns with the identity permutation
      std::vector<int> cperm(A.cols());
      int k = 0;
      for(std::vector<int>::iterator it = cperm.begin(); it != cperm.end(); ++it){
        *it = k++;
      }

      // a permutation of row j applied to all rows must not change rows 1 to j-1
      // hence, a permutation of row j can only permute entries in blocks for
      // which the columns above the entries in a block are equal
      // 
      // initially we can permute the complete row, so we have a single
      // block of size equal to the number of columns
      std::vector<int> blocks(1);
      blocks[0] = A.cols()-1;

      
      DistanceMatrixPermutation dmp(rperm,cperm,blocks,sequence(0,A.rows()));
      std::vector<DistanceMatrixPermutation> dmp_list_in;
      dmp_list_in.push_back(dmp);

      // generate all possible permutations of the distance matrix that lead to a lex max matrix.
      for ( int i = 0; i < A.rows(); ++i ) {
	std::vector<DistanceMatrixPermutation> dmp_list;
	Vector<Integer> base(A.cols());
	for ( std::vector<DistanceMatrixPermutation>::const_iterator dmp_it = dmp_list_in.begin(); dmp_it != dmp_list_in.end(); ++dmp_it ) 
	  get_all_permutations_for_row(0,*dmp_it,A,dmp_list,base);
	dmp_list_in = dmp_list;
      }


      Matrix<Integer> VR = V - repeat_row(V[0],V.rows());
      Matrix<Integer> W = common::flint::HermiteNormalForm(dmp_list_in[0].apply_vertex_permutation(VR));

      for ( int j = 1; j < W.rows(); ++j ) {
	VR = V - repeat_row(V[j],V.rows());
	for ( int i = 0; i < dmp_list_in.size(); ++i ) {
	  Matrix<Integer> U = common::flint::HermiteNormalForm(dmp_list_in[i].apply_vertex_permutation(VR));
	  if ( T(U) < T(W) ) { W = U; }
#ifdef DEBUG
	  cout << "printing permutation " << i << endl;
	  //	  cout << dmp_list_in[i] << endl;
	  //	  cout << dmp_list_in[i].apply_permutation(A) << endl;
	  cout << W << endl;
#endif
	}
      }
      
      return (ones_vector<Integer>(W.rows()))|W;
    }

    bool lattice_isomorphic ( const perl::Object & p, const perl::Object & q ) {

      Matrix<Integer> pm = affine_lattice_normal_form(p);
      Matrix<Integer> qm = affine_lattice_normal_form(q);
      
      return pm == qm;
    }
    
    
    UserFunction4perl("# @category Lattice Polytopes\n"
      "# Takes a lattice polytope and computes the affine normal form of the vertices as defined in ."
      "\tA. Grinis, Kasprzyk, Normal Forms of Polytopes, arxiv:1301.6641, Jan 2013\n"
      "# @param LatticePolytope P\n"
      "# @return Matrix<Integer>", &affine_lattice_normal_form, "affine_lattice_normal_form(Polytope)");
    
    UserFunction4perl("# @category Lattice Polytopes\n"
      "# Checks whether two lattice polytopes are lattice isomoprhic by comparing their normal forms.\n"
      "# @param LatticePolytope P\n"
      "# @param LatticePolytope Q\n"
      "# @return bool", &lattice_isomorphic, "lattice_isomorphic(Polytope, Polytope)");
  }
}
