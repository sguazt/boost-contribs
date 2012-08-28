#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cstddef>
#include <iostream>

namespace ublas = boost::numeric::ublas;

typedef int value_type;

void test_integral_vector()
{
	std::cout << "[test_integral_vector] BEGIN" << std::endl;

	std::size_t n(5);

	ublas::zero_vector<value_type> in(n);

	std::cerr << "[test_integral_vector] Input Vector: " << in << std::endl;

	// Test copy-construction
	try
	{
		ublas::vector<value_type> out(in);

		std::cout << "[test_integral_vector] Copy-construction succeeded." << std::endl;
		std::cout << "[test_integral_vector] Copy-constructed Output Vector: " << out << std::endl;
	}
	catch (...)
	{
		std::cout << "[test_integral_vector] Copy-construction failed." << std::endl;
	}

	// Test copy-assignement
	try
	{
		ublas::vector<value_type> out;
		out = in;

		std::cout << "[test_integral_vector] Copy-construction succeeded." << std::endl;
		std::cout << "[test_integral_vector] Copy-constructed Output Vector: " << out << std::endl;
	}
	catch (...)
	{
		std::cout << "[test_integral_vector] Copy-assignement failed." << std::endl;
	}

	std::cout << "[test_integral_vector] END" << std::endl;
}


void test_integral_matrix()
{
	std::cout << "[test_integral_matrix] BEGIN" << std::endl;

	std::size_t n(5);

	ublas::triangular_matrix<value_type,ublas::lower> IN(n,n);

	// Initialize the input matrix
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j <= i; ++j)
		{
			IN(i,j) = 2*i+j;
		}
	}

	std::cout << "[test_integral_matrix] Input Matrix: " << IN << std::endl;

	// Test copy-construction
	try
	{
		ublas::matrix<value_type> OUT(IN);

		std::cout << "[test_integral_matrix] Copy-construction succeeded." << std::endl;
		std::cout << "[test_integral_matrix] Copy-constructed Output Matrix: " << OUT << std::endl;
	}
	catch (...)
	{
		std::cout << "[test_integral_matrix] Copy-construction failed." << std::endl;
	}

	// Test copy-assignement
	try
	{
		ublas::matrix<value_type> OUT(n,n);
		OUT = IN;

		std::cout << "[test_integral_matrix] Copy-assignement succeeded." << std::endl;
		std::cout << "[test_integral_matrix] Copy-assigned Output Matrix: " << OUT << std::endl;
	}
	catch (...)
	{
		std::cout << "[test_integral_matrix] Copy-assignement failed." << std::endl;
	}

	std::cout << "[test_integral_matrix] END" << std::endl;

}


int main()
{
	test_integral_vector();
	test_integral_matrix();
}

