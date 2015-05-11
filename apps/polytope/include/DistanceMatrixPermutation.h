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

#include <vector>
#include <iostream>
#include <polymake/Set.h>
#include <polymake/Matrix.h>
#include <polymake/Integer.h>
#include <polymake/Rational.h>

namespace polymake {
  namespace polytope {

    // keeps a partially constructed max distance matrix
    // we actually only store
    //   the column permutation
    //   the permutation for the rows already included in the permuted distance matrix
    //   the rows we haven't looked at so far
    //   the decomposition of the column set into blocks, any further row perm must stay inside the blocks so that the already constructed rows don't change
    //     i.e. if rows 1, ..., k are already permuted to an initial max distance matrix,
    //          then any permutation we apply to the next row must leave the first k rows invariant (the entries above the row must be constant in each block)
    //
    // FIXME we should check whether it speed up computations if we store the actual permuted matrix instead of recomputing permuted rows if we need them
    //       for some polytopes we need to create quite a lot of instances of these permutations though
    class DistanceMatrixPermutation {
      
    private :
      std::vector<int> rperm;        // the permutation applied to the rows of the matrix,
                                     // maybe filled only partially if used for constructing possible permutations 
      std::vector<int> cperm;        // the permutation applied to the columns of the matrix
      std::vector<int> blocks;       // once some rows of the matrix are already permuted into lex max format
                                     // permutations of columns of further rows are only allowed among those columns that
                                     // are constant on the already permuted rows
                                     // the vector blocks lists these columns in the form
                                     // i1,i2,...,ik, where the already permuted columns are constant on 0..i1, i1+1..i2, i2+1..i3, ...
      Set<int> still_available_rows; // a list of row indices that are not yet included in the row permutation rperm
                                     // i.e. the candidate set for the next row added to a potentially max permutation

      /*
	Note that permutations are given in inverted form: cperm[i] is the index that should go into position i in the permuted vector
      */
  
    public :
      // Constructors
    DistanceMatrixPermutation(const std::vector<int> rperm_in, const std::vector<int> cperm_in, const std::vector<int> blocks_in, const Set<int> sar)
      : rperm(rperm_in), cperm(cperm_in), blocks(blocks_in), still_available_rows(sar) { }
      
    DistanceMatrixPermutation(const DistanceMatrixPermutation & dmp)
      : rperm(dmp.get_rperm()), cperm(dmp.get_cperm()), blocks(dmp.get_blocks()), still_available_rows(dmp.get_still_available_rows()) { } 
      
      
      // return private vars
      const std::vector<int> get_rperm()                 const { return rperm; }
      const std::vector<int> get_cperm()                 const { return cperm; }
      const std::vector<int> get_blocks()                const { return blocks; }
      const Set<int>         get_still_available_rows () const { return still_available_rows; }
      
      // size functions
      const int rowsize() const { return rperm.size(); }
      const int colsize() const { return cperm.size(); }
      const int blocksize() const { return blocks.size(); }

      // apply a row/column permutation to a distance matrix
      const Matrix<Integer> apply_permutation(const Matrix<Integer>& A) {

	Matrix<Integer> B(A.rows(),A.cols());
	for ( int i = 0; i < A.cols(); ++i )
	  B.col(i) = A.col(cperm[i]);

	Matrix<Integer>C(A.rows(),A.cols());
	for ( int i = 0; i < A.rows(); ++i )
	  C.row(i) = B.row(rperm[i]);
	
	return C;
      }

      // apply a column permutation to a matrix of vertices
      const Matrix<Integer> apply_vertex_permutation(const Matrix<Integer>& A) {

	Matrix<Integer> B(A.rows(),A.cols());
	for ( int i = 0; i < A.rows(); ++i )
	  B.row(i) = A.row(cperm[i]);

	return B;
      }
      

      // printing
      friend
	std::ostream& operator<< (std::ostream & os, const DistanceMatrixPermutation & dmp) {
	
	const std::vector<int> rperm = dmp.get_rperm();
	const std::vector<int> cperm = dmp.get_cperm();
	const std::vector<int> blocks = dmp.get_blocks();    
	
	os << "row permutation: ";
	for ( int j = 0; j < dmp.rowsize(); ++j ) 
	  os << rperm[j] << " ";
	os << std::endl;
	
	os << "column permutation: ";
	for ( int j = 0; j < dmp.colsize(); ++j ) 
	  os << cperm[j] << " ";
	os << std::endl;
	
	os << "blocks: ";
	for ( int j = 0; j < dmp.blocksize(); ++j ) 
	  os << blocks[j] << " ";
	os << std::endl;
	
	os << "still available rows: ";
	wrap(os) << dmp.get_still_available_rows();
	os << std::endl;
	
	return os;
      }
      
    };

  }
}
