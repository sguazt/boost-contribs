/**
 *  \file util.hpp
 *
 *  \brief Utility macros/functions for testing and debugging purpose.
 *
 *  Copyright (c) 2009, Marco Guazzone
 *
 *  Distributed under the Boost Software License, Version 1.0. (See
 *  accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 *
 * \author Marco Guazzone, marco.guazzone@gmail.com
 */

#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP


#include <iostream>


///< Expand its argument.
#define EXPAND_(x) x


///< Transform its argument into a string.
#define STRINGIFY_(x) #x


///< Concatenate its two \e string arguments.
#define JOIN_(x,y) x ## y


///< Output the message \a x if in debug-mode; otherwise output nothing.
#ifndef NDEBUG
# 	define BOOST_UBLAS_DEBUG_TRACE(x) std::cerr << "[Debug>> " << EXPAND_(x) << std::endl
#else
# 	define BOOST_UBLAS_DEBUG_TRACE(x) /**/
#endif // NDEBUG


///< Define the beginning of a test suite.
#define BOOST_UBLAS_TEST_BEGIN() 	/* [BOOST_UBLAS_TEST_BEGIN] */ \
									{ /* Begin of Test Suite */ \
										unsigned int test_fails_(0) \
									/* [/BOOST_UBLAS_TEST_BEGIN] */


///< Define a test case \a x inside the current test suite.
#define BOOST_UBLAS_TEST_DEF(x) void EXPAND_(x)(unsigned int& test_fails_)


///< Call the test case \a x.
#define BOOST_UBLAS_TEST_DO(x) 	/* [BOOST_UBLAS_TEST_DO] */ \
								try \
								{ \
									EXPAND_(x)(test_fails_); \
								} \
								catch (std::exception& e) \
								{ \
									++test_fails_; \
									BOOST_UBLAS_TEST_ERROR( e.what() ); \
								} \
								catch (...) \
								{ \
									++test_fails_; \
								} \
								/* [/BOOST_UBLAS_TEST_DO] */


///< Define the end of a test suite.
#define BOOST_UBLAS_TEST_END() 	/* [BOOST_UBLAS_TEST_END] */ \
								if (test_fails_ > 0) \
								{ \
									std::cerr << "Number of failed tests: " << test_fails_ << std::endl; \
								} \
								else \
								{ \
									std::cerr << "No failed test" << std::endl; \
								} \
								} /* End of test suite */ \
								/* [/BOOST_UBLAS_TEST_END] */


///< Check the truth of assertion \a x.
#define BOOST_UBLAS_TEST_CHECK(x) if (!(x)) { BOOST_UBLAS_TEST_ERROR( "Failed assertion: " << STRINGIFY_(x) ); ++test_fails_; }


///< Output the error message \a x.
#define BOOST_UBLAS_TEST_ERROR(x) std::cerr << "[Error>> " << EXPAND_(x) << std::endl

#endif // TEST_UTILS_HPP
