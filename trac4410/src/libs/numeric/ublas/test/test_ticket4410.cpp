#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <complex>

namespace ublas = boost::numeric::ublas;

int main()
{
	typedef double in_value_type;
	typedef std::complex<double> out_value_type;

	const std::size_t n = 4;

	// Input Matrix
	ublas::symmetric_matrix<in_value_type,ublas::lower> IN(n,n);

	// Test copy-constructor: fails to compile without patch 4410
	ublas::symmetric_matrix<out_value_type,ublas::lower> OUT1(IN);

	// Test copy-assignement: fails to compile without patch 4410
	ublas::symmetric_matrix<out_value_type,ublas::lower> OUT2;
	OUT2 = IN;

	std::cout << "Input Matrix: " << IN << std::endl;
	std::cout << "Copy-Constructed Output Matrix: " << OUT1 << std::endl;
	std::cout << "Copy-Assigned Output Matrix: " << OUT2 << std::endl;
}
