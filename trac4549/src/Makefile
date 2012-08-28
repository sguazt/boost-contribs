#.PHONY: all patch clean

test_path=libs/numeric/ublas/test

CC=$(CXX)
BOOST_ROOT=$(HOME)/projects/svn/boost-trunk
CXXFLAGS=-Wall -ansi -pedantic -O0 -g -I$(BOOST_ROOT)

all: $(test_path)/test_ticket4549

$(test_path)/test_ticket4549: $(test_path)/test_ticket4549.o

#patch:
#	cp boost/numeric/ublas/detail/config.hpp.orig boost/numeric/ublas/detail/config.hpp
#	cp boost/numeric/ublas/detail/matrix_assign.hpp.orig boost/numeric/ublas/detail/matrix_assign.hpp
#	cp boost/numeric/ublas/detail/vector_assign.hpp.orig boost/numeric/ublas/detail/vector_assign.hpp
#	patch -p2 < vector_matrix_assign-allow_integral_types.patch

clean:
	rm -f $(test_path)/test_ticket4549 $(test_path)/test_ticket4549.o
