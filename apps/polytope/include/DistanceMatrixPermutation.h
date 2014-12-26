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

class DistanceMatrixPermutation {
  
 private :
  std::vector<int> rperm;
  std::vector<int> cperm;
  std::vector<int> blocks;
  
  
 public :
  DistanceMatrixPermutation(const std::vector<int> rperm_in, const std::vector<int> cperm_in, const std::vector<int> blocks_in) : rperm(rperm_in), cperm(cperm_in), blocks(blocks_in) { }
  DistanceMatrixPermutation(const DistanceMatrixPermutation & dmp) : rperm(dmp.get_rperm()), cperm(dmp.get_cperm()), blocks(dmp.get_blocks()) { } 
  
  const std::vector<int> get_rperm()  const { return rperm; }
  const std::vector<int> get_cperm()  const { return cperm; }
  const std::vector<int> get_blocks() const { return blocks; }


  friend
    std::ostream& operator<< (std::ostream & os, const DistanceMatrixPermutation & dpm) {

    os << "row permutation: ";
    for ( int j = 0; j < dpm.get_rperm().size(); ++j ) 
      os << dpm.get_rperm()[j] << " ";
    os << std::endl;

    os << "column permutation: ";
    for ( int j = 0; j < dpm.get_cperm().size(); ++j ) 
      os << dpm.get_cperm()[j] << " ";
    os << std::endl;

    os << "blocks: ";
    for ( int j = 0; j < dpm.get_blocks().size(); ++j ) 
      os << dpm.get_blocks()[j] << " ";
    os << std::endl;

    return os;
   }

};
