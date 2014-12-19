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


namespace polymake { 

  namespace polytope {

    
    Matrix<Integer> lattice_normal_form ( const Matrix<Integer> & A ) {
      
      return common::HermiteNormalForm(A);
    } 
    
    UserFunction4perl(" ", &lattice_normal_form, "lattice_normal_form($)");
        
  }
}
