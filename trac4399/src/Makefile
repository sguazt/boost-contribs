src_path=.
test_path=libs/numeric/ublas/test
apidoc_path=libs/numeric/ublas/doc/api

CXXFLAGS=-Wall -Wextra -pedantic -ansi -I$(src_path) -I$(HOME)/projects/svn/boost-trunk

LDFLAGS=-lm

CC=$(CXX)
CLEANER=rm -rf
DOXYGEN=doxygen


all: 	$(test_path)/size

$(test_path)/size: $(test_path)/size.o

apidoc:
	$(DOXYGEN) Doxyfile

clean:
	$(CLEANER)	$(test_path)/size $(test_path)/size.o \
				$(apidoc_path)

