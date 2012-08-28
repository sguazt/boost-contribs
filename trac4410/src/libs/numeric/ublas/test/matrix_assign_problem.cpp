#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace ublas = boost::numeric::ublas;

int main()
{
    typedef double in_value_type;
    typedef std::complex<double> out_value_type;

    const std::size_t n = 4;

    ublas::matrix<in_value_type> IN(n,n);
    ublas::matrix<out_value_type> OUT(IN); // COMPILE
    OUT = IN; // COMPILE

    ublas::symmetric_matrix<in_value_type,ublas::lower> sym_IN(n,n);
    ublas::symmetric_matrix<out_value_type,ublas::lower> sym_OUT(sym_IN); // NOT COMPILE
    OUT = IN; // NOT COMPILE
}
