/**
 *  \file generalied_diagonal_matrix.cpp
 *
 *  \brief Test suite for the \c generalied_diagonal_matrix matrix container.
 *
 *  Copyright (c) 2009, Marco Guazzone
 *
 *  Distributed under the Boost Software License, Version 1.0. (See
 *  accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 *
 *  \author Marco Guazzone, marco.guazzone@gmail.com
 */

#include <boost/numeric/ublas/container/generalized_diagonal_matrix.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <iostream>
#include "libs/numeric/ublas/test/utils.hpp"


static const double TOL(1.0e-5); ///< Tolerance for real numbers comparison.


//@{ Construction //////////////////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;

	BOOST_UBLAS_DEBUG_TRACE( "A(0,0) " << A(0,0) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,0) - 0.555950) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,1) " << A(1,1) << " ==> " << 0.830123 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,1) - 0.830123) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,2) " << A(2,2) << " ==> " << 0.216504 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,2) - 0.216504) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,3) " << A(3,3) << " ==> " << 0.450332 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,3) - 0.450332) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,1) " << A(0,1) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,1) - 0.274690) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,2) " << A(1,2) << " ==> " << 0.891726 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,2) - 0.891726) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,3) " << A(2,3) << " ==> " << 0.883152 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,3) - 0.883152) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,2) " << A(0,2) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,2) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,3) " << A(1,3) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,3) - 0.895283) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,3) " << A(0,3) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,3) - 0.798938) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(1,0) " << A(1,0) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,0) - 0.108929) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,1) " << A(2,1) << " ==> " << 0.973234 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,1) - 0.973234) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,2) " << A(3,2) << " ==> " << 0.231751 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,2) - 0.231751) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(2,0) " << A(2,0) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,0) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,1) " << A(3,1) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,1) - 0.675382) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(3,0) " << A(3,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,0) - 0.023787) <= TOL );
}


//@} Construction //////////////////////////////////////////////////////////////

//@{ Column-major construction /////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;

	BOOST_UBLAS_DEBUG_TRACE( "A(0,0) " << A(0,0) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,0) - 0.555950) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,1) " << A(1,1) << " ==> " << 0.830123 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,1) - 0.830123) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,2) " << A(2,2) << " ==> " << 0.216504 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,2) - 0.216504) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,3) " << A(3,3) << " ==> " << 0.450332 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,3) - 0.450332) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, 1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,1) " << A(0,1) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,1) - 0.274690) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,2) " << A(1,2) << " ==> " << 0.891726 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,2) - 0.891726) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,3) " << A(2,3) << " ==> " << 0.883152 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,3) - 0.883152) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, 2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,2) " << A(0,2) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,2) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,3) " << A(1,3) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,3) - 0.895283) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, 3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,3) " << A(0,3) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,3) - 0.798938) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, -1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(1,0) " << A(1,0) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,0) - 0.108929) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,1) " << A(2,1) << " ==> " << 0.973234 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,1) - 0.973234) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,2) " << A(3,2) << " ==> " << 0.231751 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,2) - 0.231751) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(2,0) " << A(2,0) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,0) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,1) " << A(3,1) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,1) - 0.675382) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_col_major )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Column Major" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;

	matrix_type A(4, -3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(3,0) " << A(3,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,0) - 0.023787) <= TOL );
}


//@} Column-major construction /////////////////////////////////////////////////

//@{ Rectangular Construction //////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_hrect_main_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,0);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332; /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,0) " << A(0,0) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,0) - 0.555950) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,1) " << A(1,1) << " ==> " << 0.830123 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,1) - 0.830123) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,2) " << A(2,2) << " ==> " << 0.216504 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,2) - 0.216504) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,3) " << A(3,3) << " ==> " << 0.450332 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,3) - 0.450332) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            A(3,4) = 0.555950; /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,1) " << A(0,1) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,1) - 0.274690) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,2) " << A(1,2) << " ==> " << 0.891726 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,2) - 0.891726) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,3) " << A(2,3) << " ==> " << 0.883152 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,3) - 0.883152) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,4) " << A(3,4) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,4) - 0.555950) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            A(2,4) = 0.555950; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(3,5) = 0.274690; /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,2) " << A(0,2) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,2) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,3) " << A(1,3) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,3) - 0.895283) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,4) " << A(2,4) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,4) - 0.555950) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,5) " << A(3,5) << " ==> " << 0.274690 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,5) - 0.274690) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            A(1,4) = 0.540605; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(2,5) = 0.895283; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(3,6) = 0.555950;

	BOOST_UBLAS_DEBUG_TRACE( "A(0,3) " << A(0,3) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,3) - 0.798938) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,4) " << A(1,4) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,4) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,5) " << A(2,5) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,5) - 0.895283) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,6) " << A(3,6) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,6) - 0.555950) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up4_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Fourth Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,4);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            A(0,4) = 0.798938; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(1,5) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(2,6) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,4) " << A(0,4) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,4) - 0.798938) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,5) " << A(1,5) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,5) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,6) " << A(2,6) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,6) - 0.895283) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up5_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Fifth Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,5);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(0,5) = 0.798938; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(1,6) = 0.540605;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,5) " << A(0,5) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,5) - 0.798938) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,6) " << A(1,6) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,6) - 0.540605) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_up6_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Sixth Upper Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,6);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            A(0,6) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,6) " << A(0,6) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,6) - 0.798938) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_low1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,-1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(1,0) " << A(1,0) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,0) - 0.108929) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,1) " << A(2,1) << " ==> " << 0.973234 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,1) - 0.973234) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,2) " << A(3,2) << " ==> " << 0.231751 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,2) - 0.231751) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_low2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,-2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(2,0) " << A(2,0) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,0) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,1) " << A(3,1) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,1) - 0.675382) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_hrect_low3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Horizontal Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4,7,-3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(3,0) " << A(3,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,0) - 0.023787) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_main_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,0);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,0) " << A(0,0) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,0) - 0.555950) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,1) " << A(1,1) << " ==> " << 0.830123 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,1) - 0.830123) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,2) " << A(2,2) << " ==> " << 0.216504 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,2) - 0.216504) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,3) " << A(3,3) << " ==> " << 0.450332 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,3) - 0.450332) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_up1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,1) " << A(0,1) << " ==> " << 0.555950 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,1) - 0.274690) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,2) " << A(1,2) << " ==> " << 0.891726 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,2) - 0.891726) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,3) " << A(2,3) << " ==> " << 0.883152 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,3) - 0.883152) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_up2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,2) " << A(0,2) << " ==> " << 0.540605 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,2) - 0.540605) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(1,3) " << A(1,3) << " ==> " << 0.895283 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,3) - 0.895283) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_up3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(0,3) " << A(0,3) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(0,3) - 0.798938) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low1_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(4,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(1,0) " << A(1,0) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(1,0) - 0.108929) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(2,1) " << A(2,1) << " ==> " << 0.973234 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,1) - 0.973234) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,2) " << A(3,2) << " ==> " << 0.231751 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,2) - 0.231751) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(4,3) " << A(4,3) << " ==> " << 0.798938 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(4,3) - 0.798938) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low2_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 0.108929; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(5,3) = 0.973234;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(2,0) " << A(2,0) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(2,0) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(3,1) " << A(3,1) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,1) - 0.675382) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(4,2) " << A(4,2) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(4,2) - 0.108929) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(5,3) " << A(5,3) << " ==> " << 0.973234 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(5,3) - 0.973234) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low3_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(4,1) = 0.948014; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(5,2) = 0.675382; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(6,3) = 0.108929;

	BOOST_UBLAS_DEBUG_TRACE( "A(3,0) " << A(3,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(3,0) - 0.023787) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(4,1) " << A(4,1) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(4,1) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(5,2) " << A(5,2) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(5,2) - 0.675382) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(6,3) " << A(6,3) << " ==> " << 0.108929 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(6,3) - 0.108929) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low4_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Fourth Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-4);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(4,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(5,1) = 0.948014; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(6,2) = 0.675382; /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(4,0) " << A(4,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(4,0) - 0.023787) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(5,1) " << A(5,1) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(5,1) - 0.948014) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(6,2) " << A(6,2) << " ==> " << 0.675382 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(6,2) - 0.675382) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low5_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Fifth Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-5);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(5,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(6,1) = 0.948014; /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(5,0) " << A(5,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(5,0) - 0.023787) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "A(6,1) " << A(6,1) << " ==> " << 0.948014 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(6,1) - 0.948014) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_vrect_low6_diagonal )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Sixth Lower Diagonal -- Vertical Rectangular Matrix" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(7,4,-6);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(6,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	BOOST_UBLAS_DEBUG_TRACE( "A(6,0) " << A(6,0) << " ==> " << 0.023787 );
	BOOST_UBLAS_TEST_CHECK( std::fabs(A(6,0) - 0.023787) <= TOL );
}


//@} Rectangular Construction //////////////////////////////////////////////////

//@{ Row-by-Column Iteration ///////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_row_col_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Row-Col Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator1 row_it = A.begin1();
		row_it != A.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *col_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - A(row,col)) <= TOL );
		}
	}
}


//@} Row-by-Column Iteration ///////////////////////////////////////////////////

//@{ Column-by-Row Iteration ///////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_col_row_iteration )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Col-Row Iteration" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	for (
		matrix_type::const_iterator2 col_it = A.begin2();
		col_it != A.end2();
		++col_it
	) {
		for (
			matrix_type::const_iterator1 row_it = col_it.begin();
			row_it != col_it.end();
			++row_it
		) {
			matrix_type::size_type row(row_it.index1());
			matrix_type::size_type col(row_it.index2());

			BOOST_UBLAS_DEBUG_TRACE( "A(" << row << "," << col << ") " << *row_it << " ==> " << A(row,col) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*row_it - A(row,col)) <= TOL );
		}
	}
}


//@} Column-by-Row Iteration ///////////////////////////////////////////////////

//@{ Copy-Construction /////////////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4);

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,0) " << B(0,0) << " ==> " << A(0,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,0) - A(0,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,1) " << B(1,1) << " ==> " << A(1,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,1) - A(1,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,2) " << B(2,2) << " ==> " << A(2,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,2) - A(2,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,3) " << B(3,3) << " ==> " << A(3,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,3) - A(3,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 1);

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,1) " << B(0,1) << " ==> " << A(0,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,1) - A(0,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,2) " << B(1,2) << " ==> " << A(1,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,2) - A(1,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,3) " << B(2,3) << " ==> " << A(2,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,3) - A(2,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 2);

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,2) " << B(0,2) << " ==> " << A(0,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,2) - A(0,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,3) " << B(1,3) << " ==> " << A(1,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,3) - A(1,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, 3);

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,3) " << B(0,3) << " ==> " << A(0,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,3) - A(0,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -1);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(1,0) " << B(1,0) << " ==> " << A(1,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,0) - A(1,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,1) " << B(2,1) << " ==> " << A(2,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,1) - A(2,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,2) " << B(3,2) << " ==> " << A(3,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,2) - A(3,2)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(2,0) " << B(2,0) << " ==> " << A(2,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,0) - A(2,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,1) " << B(3,1) << " ==> " << A(3,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,1) - A(3,1)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(4, -3);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */

	matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(3,0) " << B(3,0) << " ==> " << A(3,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,0) - A(3,0)) <= TOL );
}


//@} Copy-Construction /////////////////////////////////////////////////////////

//@{ Matrix-Copy-Construction //////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	A(0,0) = 0.555950; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(1,1) = 0.830123; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(2,2) = 0.216504; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(3,3) = 0.450332;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	result_matrix_type B(A);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,0) " << B(0,0) << " ==> " << A(0,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,0) - A(0,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,1) " << B(1,1) << " ==> " << A(1,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,1) - A(1,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,2) " << B(2,2) << " ==> " << A(2,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,2) - A(2,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,3) " << B(3,3) << " ==> " << A(3,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,3) - A(3,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            A(0,1) = 0.274690; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(1,2) = 0.891726; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(2,3) = 0.883152;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	result_matrix_type B(A, 1);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,1) " << B(0,1) << " ==> " << A(0,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,1) - A(0,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,2) " << B(1,2) << " ==> " << A(1,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,2) - A(1,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,3) " << B(2,3) << " ==> " << A(2,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,3) - A(2,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            A(0,2) = 0.540605; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(1,3) = 0.895283;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	result_matrix_type B(A, 2);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,2) " << B(0,2) << " ==> " << A(0,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,2) - A(0,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,3) " << B(1,3) << " ==> " << A(1,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,3) - A(1,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            /* 0 */            A(0,3) = 0.798938;
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */

	result_matrix_type B(A, 3);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,3) " << B(0,3) << " ==> " << A(0,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,3) - A(0,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(1,0) = 0.108929; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(2,1) = 0.973234; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(3,2) = 0.231751; /* 0 */
	/* 0 */            /* 0 */            /* 0 */            A(4,3) = 1.450332;

	result_matrix_type B(A, -1);

	BOOST_UBLAS_DEBUG_TRACE( "B(1,0) " << B(1,0) << " ==> " << A(1,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,0) - A(1,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,1) " << B(2,1) << " ==> " << A(2,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,1) - A(2,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,2) " << B(3,2) << " ==> " << A(3,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,2) - A(3,2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(4,3) " << B(4,3) << " ==> " << A(4,3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(4,3) - A(4,3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	result_matrix_type B(A, -2);

	BOOST_UBLAS_DEBUG_TRACE( "B(2,0) " << B(2,0) << " ==> " << A(2,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,0) - A(2,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,1) " << B(3,1) << " ==> " << A(3,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,1) - A(3,1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(4,2) " << B(4,2) << " ==> " << A(4,2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(4,2) - A(4,2)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(3,0) = 0.023787; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(4,1) = 1.675382; /* 0 */            /* 0 */

	result_matrix_type B(A, -3);

	BOOST_UBLAS_DEBUG_TRACE( "B(3,0) " << B(3,0) << " ==> " << A(3,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,0) - A(3,0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(4,1) " << B(4,1) << " ==> " << A(4,1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(4,1) - A(4,1)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low4_diagonal_mat_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Fourth Lower Diagonal -- Matrix-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::matrix<value_type> source_matrix_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> result_matrix_type;

	source_matrix_type A(5, 4, value_type(0));

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(4,0) = 1.023787; /* 0 */            /* 0 */            /* 0 */

	result_matrix_type B(A, -4);

	BOOST_UBLAS_DEBUG_TRACE( "B(4,0) " << B(4,0) << " ==> " << A(4,0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(4,0) - A(4,0)) <= TOL );
}


//@} Matrix-Copy-Construction //////////////////////////////////////////////////

//@{ Vector-Copy-Construction //////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_main_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Main Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(4);

	v(0) = 0.555950;
	v(1) = 0.830123;
	v(2) = 0.216504;
	v(3) = 0.450332;

	matrix_type B(v);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,0) " << B(0,0) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,0) - v(0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,1) " << B(1,1) << " ==> " << v(1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,1) - v(1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,2) " << B(2,2) << " ==> " << v(2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,2) - v(2)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,3) " << B(3,3) << " ==> " << v(3) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,3) - v(3)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up1_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Upper Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(3);

	v(0) = 0.274690;
	v(1) = 0.891726;
	v(2) = 0.883152;

	matrix_type B(v, 1);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,1) " << B(0,1) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,1) - v(0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,2) " << B(1,2) << " ==> " << v(1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,2) - v(1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,3) " << B(2,3) << " ==> " << v(2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,3) - v(2)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up2_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Upper Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(2);

	v(0) = 0.540605;
	v(1) = 0.895283;

	matrix_type B(v, 2);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,2) " << B(0,2) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,2) - v(0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(1,3) " << B(1,3) << " ==> " << v(1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,3) - v(1)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_up3_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Upper Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(3);

	v(0) = 0.798938;

	matrix_type B(v, 3);

	BOOST_UBLAS_DEBUG_TRACE( "B(0,3) " << B(0,3) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(0,3) - v(0)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low1_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST First Lower Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(3);

	v(0) = 0.108929;
	v(1) = 0.973234;
	v(2) = 0.231751;

	matrix_type B(v, -1);

	BOOST_UBLAS_DEBUG_TRACE( "B(1,0) " << B(1,0) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(1,0) - v(0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(2,1) " << B(2,1) << " ==> " << v(1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,1) - v(1)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,2) " << B(3,2) << " ==> " << v(2) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,2) - v(2)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low2_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Second Lower Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(2);

	v(0) = 0.948014;
	v(1) = 0.675382;

	matrix_type B(v, -2);

	BOOST_UBLAS_DEBUG_TRACE( "B(2,0) " << B(2,0) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(2,0) - v(0)) <= TOL );
	BOOST_UBLAS_DEBUG_TRACE( "B(3,1) " << B(3,1) << " ==> " << v(1) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,1) - v(1)) <= TOL );
}


BOOST_UBLAS_TEST_DEF( test_low3_diagonal_vec_copy )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Third Lower Diagonal -- Vector-Copy-Construction" );

	typedef double value_type;
	typedef boost::numeric::ublas::vector<value_type> vector_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	vector_type v(1);

	v(0) = 0.023787;

	matrix_type B(v, -3);

	BOOST_UBLAS_DEBUG_TRACE( "B(3,0) " << B(3,0) << " ==> " << v(0) );
	BOOST_UBLAS_TEST_CHECK( std::fabs(B(3,0) - v(0)) <= TOL );
}


//@} Vector-Copy-Construction //////////////////////////////////////////////////

//@{ Matrix Operations /////////////////////////////////////////////////////////


BOOST_UBLAS_TEST_DEF( test_op_transpose )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Transpose" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type C(4, 5, 2);
	C = boost::numeric::ublas::trans(A);

	for (
		matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			const matrix_type::value_type a = static_cast<matrix_type const&>(A)(col,row); // FIXME: This type of cast is very boring!!

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << a );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - a) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_op_sum_dense )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Generalized Diagonal + Dense" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type1;
	typedef boost::numeric::ublas::matrix<value_type> matrix_type2;
	typedef boost::numeric::ublas::matrix<value_type> result_matrix_type;

	matrix_type1 A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type2 B(5,4);

	B(0,0) = 0.555950; B(0,1) = 0.274690; B(0,2) = 0.540605; B(0,3) = 0.798938;
	B(1,0) = 0.108929; B(1,1) = 0.830123; B(1,2) = 0.891726; B(1,3) = 0.895283;
	B(2,0) = 0.948014; B(2,1) = 0.973234; B(2,2) = 0.216504; B(2,3) = 0.883152;
	B(3,0) = 0.023787; B(3,1) = 0.675382; B(3,2) = 0.231751; B(3,3) = 0.450332;
	B(4,0) = 1.023787; B(4,1) = 1.675382; B(4,2) = 1.231751; B(4,3) = 1.450332;

	result_matrix_type C = A + B;

	for (
		result_matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			result_matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			result_matrix_type::size_type row(col_it.index1());
			result_matrix_type::size_type col(col_it.index2());

			const matrix_type1::value_type a = static_cast<matrix_type1 const&>(A)(row,col); // FIXME: This type of cast is very boring!!
			const matrix_type2::value_type b = B(row,col);

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << (a+b) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - (a+b)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_op_diff_dense )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Generalized Diagonal - Dense" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type1;
	typedef boost::numeric::ublas::matrix<value_type> matrix_type2;
	typedef boost::numeric::ublas::matrix<value_type> result_matrix_type;

	matrix_type1 A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type2 B(5,4);

	B(0,0) = 0.555950; B(0,1) = 0.274690; B(0,2) = 0.540605; B(0,3) = 0.798938;
	B(1,0) = 0.108929; B(1,1) = 0.830123; B(1,2) = 0.891726; B(1,3) = 0.895283;
	B(2,0) = 0.948014; B(2,1) = 0.973234; B(2,2) = 0.216504; B(2,3) = 0.883152;
	B(3,0) = 0.023787; B(3,1) = 0.675382; B(3,2) = 0.231751; B(3,3) = 0.450332;
	B(4,0) = 1.023787; B(4,1) = 1.675382; B(4,2) = 1.231751; B(4,3) = 1.450332;

	result_matrix_type C = A - B;

	for (
		result_matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			result_matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			result_matrix_type::size_type row(col_it.index1());
			result_matrix_type::size_type col(col_it.index2());

			const matrix_type1::value_type a = static_cast<matrix_type1 const&>(A)(row,col); // FIXME: This type of cast is very boring!!
			const matrix_type2::value_type b = B(row,col);

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << (a-b) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - (a-b)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_op_prod )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Generalized Diagonal * Generalized Diagonal" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type;

	matrix_type A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type B(4, 3, 1);

	/* 0 */            B(0,1) = 0.274690; /* 0 */
	/* 0 */            /* 0 */            B(1,2) = 0.891726;
	/* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */

	matrix_type T(5, 3, -1);

	/* 0 */            /* 0 */            /* 0 */
	T(1,0) = 0;        /* 0 */            /* 0 */
	/* 0 */            T(2,1) = 0.260410; /* 0 */
	/* 0 */            /* 0 */            T(3,2) = 0.602256;
	/* 0 */            /* 0 */            /* 0 */

	matrix_type C(5, 3, -1);
	C = boost::numeric::ublas::prod(A, B);

	for (
		matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			matrix_type::size_type row(col_it.index1());
			matrix_type::size_type col(col_it.index2());

			const matrix_type::value_type t = static_cast<matrix_type const&>(T)(row,col); // FIXME: This type of cast is very boring!!

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << t );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - t) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_op_element_prod_dense )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Generalized Diagonal .* Dense" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type1;
	typedef boost::numeric::ublas::matrix<value_type> matrix_type2;
	typedef boost::numeric::ublas::matrix<value_type> result_matrix_type;

	matrix_type1 A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type2 B(5,4);

	B(0,0) = 0.555950; B(0,1) = 0.274690; B(0,2) = 0.540605; B(0,3) = 0.798938;
	B(1,0) = 0.108929; B(1,1) = 0.830123; B(1,2) = 0.891726; B(1,3) = 0.895283;
	B(2,0) = 0.948014; B(2,1) = 0.973234; B(2,2) = 0.216504; B(2,3) = 0.883152;
	B(3,0) = 0.023787; B(3,1) = 0.675382; B(3,2) = 0.231751; B(3,3) = 0.450332;
	B(4,0) = 1.023787; B(4,1) = 1.675382; B(4,2) = 1.231751; B(4,3) = 1.450332;

	result_matrix_type C = boost::numeric::ublas::element_prod(A, B);

	for (
		result_matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			result_matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			result_matrix_type::size_type row(col_it.index1());
			result_matrix_type::size_type col(col_it.index2());

			const matrix_type1::value_type a = static_cast<matrix_type1 const&>(A)(row,col); // FIXME: This type of cast is very boring!!
			const matrix_type2::value_type b = B(row,col);

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << (a*b) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - (a*b)) <= TOL );
		}
	}
}


BOOST_UBLAS_TEST_DEF( test_op_element_div_dense )
{
	BOOST_UBLAS_DEBUG_TRACE( "TEST Operations -- Generalized Diagonal ./ Dense" );

	typedef double value_type;
	typedef boost::numeric::ublas::generalized_diagonal_matrix<value_type> matrix_type1;
	typedef boost::numeric::ublas::matrix<value_type> matrix_type2;
	typedef boost::numeric::ublas::matrix<value_type> result_matrix_type;

	matrix_type1 A(5, 4, -2);

	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	/* 0 */            /* 0 */            /* 0 */            /* 0 */
	A(2,0) = 0.948014; /* 0 */            /* 0 */            /* 0 */
	/* 0 */            A(3,1) = 0.675382; /* 0 */            /* 0 */
	/* 0 */            /* 0 */            A(4,2) = 1.231751; /* 0 */

	matrix_type2 B(5,4);

	B(0,0) = 0.555950; B(0,1) = 0.274690; B(0,2) = 0.540605; B(0,3) = 0.798938;
	B(1,0) = 0.108929; B(1,1) = 0.830123; B(1,2) = 0.891726; B(1,3) = 0.895283;
	B(2,0) = 0.948014; B(2,1) = 0.973234; B(2,2) = 0.216504; B(2,3) = 0.883152;
	B(3,0) = 0.023787; B(3,1) = 0.675382; B(3,2) = 0.231751; B(3,3) = 0.450332;
	B(4,0) = 1.023787; B(4,1) = 1.675382; B(4,2) = 1.231751; B(4,3) = 1.450332;

	result_matrix_type C = boost::numeric::ublas::element_div(A, B);

	for (
		result_matrix_type::const_iterator1 row_it = C.begin1();
		row_it != C.end1();
		++row_it
	) {
		for (
			result_matrix_type::const_iterator2 col_it = row_it.begin();
			col_it != row_it.end();
			++col_it
		) {
			result_matrix_type::size_type row(col_it.index1());
			result_matrix_type::size_type col(col_it.index2());

			const matrix_type1::value_type a = static_cast<matrix_type1 const&>(A)(row,col); // FIXME: This type of cast is very boring!!
			const matrix_type2::value_type b = B(row,col);

			BOOST_UBLAS_DEBUG_TRACE( "C(" << row << "," << col << ") " << *col_it << " ==> " << (a/b) );
			BOOST_UBLAS_TEST_CHECK( std::fabs(*col_it - (a/b)) <= TOL );
		}
	}
}


//@} Matrix Operations /////////////////////////////////////////////////////////


int main()
{
	BOOST_UBLAS_TEST_BEGIN();

	// Matrix construction tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal );

	// Column-major matrix construction tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_col_major );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_col_major );

	// Rectangular matrix construction tests
	BOOST_UBLAS_TEST_DO( test_hrect_main_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up1_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up2_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up3_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up4_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up5_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_up6_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_low1_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_low2_diagonal );
	BOOST_UBLAS_TEST_DO( test_hrect_low3_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_main_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_up1_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_up2_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_up3_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low1_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low2_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low3_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low4_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low5_diagonal );
	BOOST_UBLAS_TEST_DO( test_vrect_low6_diagonal );

	// Matrix row-by-column iteration tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_row_col_iteration );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_row_col_iteration );

	// Matrix column-by-row iteration tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_col_row_iteration );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_col_row_iteration );

	// Matrix copy-construction tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_copy );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_copy );

	// Matrix matrix-copy-construction tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_mat_copy );
	BOOST_UBLAS_TEST_DO( test_low4_diagonal_mat_copy );

	// Matrix vector-copy-construction tests
	BOOST_UBLAS_TEST_DO( test_main_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_up1_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_up2_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_up3_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_low1_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_low2_diagonal_vec_copy );
	BOOST_UBLAS_TEST_DO( test_low3_diagonal_vec_copy );

	// Matrix operations
	BOOST_UBLAS_TEST_DO( test_op_transpose );
	BOOST_UBLAS_TEST_DO( test_op_sum_dense );
	BOOST_UBLAS_TEST_DO( test_op_diff_dense );
	BOOST_UBLAS_TEST_DO( test_op_prod );
	BOOST_UBLAS_TEST_DO( test_op_element_prod_dense );
	BOOST_UBLAS_TEST_DO( test_op_element_div_dense );

	BOOST_UBLAS_TEST_END();
}
