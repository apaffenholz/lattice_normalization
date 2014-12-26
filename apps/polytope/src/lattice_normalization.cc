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
#include <polymake/common/flint_functions.h>
#include <polymake/polytope/DistanceMatrixPermutation.h>


namespace polymake { 

  namespace polytope {

    // returns a permutation pi for the vector v such that pi(v) 
    // is lex maixmal among all permutations of the entries of v
    std::pair<Vector<Integer>, Vector<int> > FindMaxPermutation ( const Vector<Integer> & v ) {

      const int dim = v.size();
      Vector<Integer> base(v);
      Vector<int> revperm(sequence(0,dim));

      for ( int i = 1; i < dim; ++i ) {
	int j = i;
	Integer h = base[i];
	int r = revperm[i];
	while ( j > 0 && base[j-1] < h ) {
	  base[j] = base[j-1];
	  revperm[j] = revperm[j-1];
	  j--;
	}
	base[j] = h;
	revperm[j] =r;
      }

      Vector<int> perm(dim);
      for ( int i = 0; i < dim; i++ ) {
	perm[revperm[i]] = i;
      }

      std::pair<Vector<Integer>, Vector<int> > ret;
      ret.first = base;
      ret.second = perm;

      return ret;
    }

    
    Matrix<Integer> lattice_normal_form ( const perl::Object & p ) {
      
      if ( ! p.give("LATTICE") ) 
	throw std::runtime_error("the given polytope is not a lattice polytope");

      Matrix<Integer> A = p.give("FACET_VERTEX_LATTICE_DISTANCES");

      cout << FindMaxPermutation(A.row(0)) << endl;

      return common::HermiteNormalForm(A);
    } 
    
    UserFunction4perl(" ", &lattice_normal_form, "lattice_normal_form(Polytope)");
        
  }
}
