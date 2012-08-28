src_path=.
test_path=libs/numeric/ublas/test
apidoc_path=libs/numeric/ublas/doc/api

CXXFLAGS=-Wall -Wextra -pedantic -ansi -I$(src_path)
LDFLAGS=-lm

CC=$(CXX)
CLEANER=rm -rf
DOXYGEN=doxygen


all: 	$(test_path)/diag \
		$(test_path)/generalized_diagonal_matrix

$(test_path)/diag: $(test_path)/diag.o

$(test_path)/generalized_diagonal_matrix: $(test_path)/generalized_diagonal_matrix.o

apidoc:
	mkdir -p $(apidoc_path)
	$(DOXYGEN) Doxyfile

clean:
	$(CLEANER)	$(test_path)/diag $(test_path)/diag.o \
				$(test_path)/generalized_diagonal_matrix $(test_path)/generalized_diagonal_matrix.o \
				$(apidoc_path)

