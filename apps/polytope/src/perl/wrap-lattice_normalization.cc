/* Copyright (c) 1997-2014
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

namespace polymake { namespace polytope { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer> const&) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Integer > > >() );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer> const&) );

   FunctionWrapper4perl( pm::Matrix<pm::Integer> (perl::Object const&) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0 );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Integer> (perl::Object const&) );

   FunctionWrapper4perl( bool (perl::Object const&, perl::Object const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1 );
   }
   FunctionWrapperInstance4perl( bool (perl::Object const&, perl::Object const&) );

   FunctionWrapper4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer>&, pm::Matrix<pm::Integer>&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< Matrix< Integer > > >(), arg1.get< perl::TryCanned< Matrix< Integer > > >() );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer>&, pm::Matrix<pm::Integer>&) );

   FunctionWrapper4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer>&, pm::Matrix<pm::Rational>&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< Matrix< Integer > > >(), arg1.get< perl::TryCanned< Matrix< Rational > > >() );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer>&, pm::Matrix<pm::Rational>&) );

   FunctionWrapper4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer> const&, pm::Matrix<pm::Integer> const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Integer > > >(), arg1.get< perl::TryCanned< const Matrix< Integer > > >() );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Integer> (pm::Matrix<pm::Integer> const&, pm::Matrix<pm::Integer> const&) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
