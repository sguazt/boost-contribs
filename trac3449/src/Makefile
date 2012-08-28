src_path=.
test_path=libs/numeric/ublas/test
apidoc_path=libs/numeric/ublas/doc/api

# For Boost 1.40
#CXXFLAGS=-Wall -Wextra -pedantic -ansi -I$(HOME)/Sys/include -I$(src_path)
# For Boost 1.37
CXXFLAGS=-Wall -Wextra -pedantic -ansi -I$(src_path)

# For Boost 1.40
#LDFLAGS=-lm -L$(HOME)/Sys/lib
# For Boost 1.37
LDFLAGS=-lm

CC=$(CXX)
CLEANER=rm -rf
DOXYGEN=doxygen


all: 	$(test_path)/begin_end \
		$(test_path)/num_columns \
		$(test_path)/num_rows \
		$(test_path)/size

$(test_path)/begin_end: $(test_path)/begin_end.o

$(test_path)/num_columns: $(test_path)/num_columns.o

$(test_path)/num_rows: $(test_path)/num_rows.o

$(test_path)/size: $(test_path)/size.o

apidoc:
	$(DOXYGEN) Doxyfile

clean:
	$(CLEANER)	$(test_path)/begin_end $(test_path)/begin_end.o \
				$(test_path)/num_rows $(test_path)/num_columns.o \
				$(test_path)/num_columns $(test_path)/num_rows.o \
				$(test_path)/size $(test_path)/size.o \
				$(apidoc_path)

