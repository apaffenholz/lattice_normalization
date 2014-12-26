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

class DistanceMatrixPermutation {
  
 private :
  std::vector<int> rperm;
  std::vector<int> cperm;
  std::vector<int> blocks;
  
  
 public :
  DistanceMatrixPermutation() { }
  DistanceMatrixPermutation(const DistanceMatrixPermutation & dmp) : rperm(dmp.get_rperm()), cperm(dmp.get_cperm()), blocks(dmp.get_blocks()) { } 
  
  const std::vector<int> get_rperm()  const { return rperm; }
  const std::vector<int> get_cperm()  const { return cperm; }
  const std::vector<int> get_blocks() const { return blocks; }

};
