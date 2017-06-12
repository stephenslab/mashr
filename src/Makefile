INCLUDE_FLAG = -I include -I . -I`Rscript -e 'cat(find.package("Rcpp"))'`/include -I`Rscript -e 'cat(find.package("RcppArmadillo"))'`/include -I`Rscript -e 'cat(R.home("include"))'`
LINK_FLAGS = -lopenblas -llapack -lm -lstdc++ # -lgomp
OPT = -O4 -DARMA_64BIT_WORD -DARMA_USE_CXX11
CXX = g++
CXXFLAGS = $(INCLUDE_FLAG) $(OPT) -std=c++11 # -fopenmp
RLIB_FLAGS = -fPIC -shared
SRC_FILES = $(filter-out test.cpp, $(wildcard *.cpp))

mashr.so: obj
	R CMD SHLIB -o mashr.so $(addsuffix .o, $(notdir $(basename $(SRC_FILES)))) $(LINK_FLAGS)

obj:
	$(CXX) $(CXXFLAGS) $(RLIB_FLAGS) -c $(SRC_FILES)

.PHONY: clean

clean:
	rm -f *.so *.o
