.PHONY: all patch clean

test_dir=libs/numeric/ublas/test
BOOST_ROOT=$(HOME)/projects/svn/boost-trunk
CXXFLAGS=-Wall -ansi -pedantic -O0 -ggdb -I$(BOOST_ROOT)

all: $(test_dir)/test_ticket4410

patch:
	cp boost/numeric/ublas/detail/matrix_assign.hpp.orig boost/numeric/ublas/detail/matrix_assign.hpp
	patch -p2 < ../patches/matrix_assign-packed_sparse_storage-type_conversion.patch

clean:
	rm -f $(test_dir)/test_ticket4410 $(test_dir)/test_ticket4410.o

realclean: clean
	rm -f boost/numeric/ublas/detail/matrix_assign.hpp
