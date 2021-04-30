# Paths to libraries
LHAPDF := /home/mbeekveld/LHAPDF
CUBA := /home/mbeekveld/Cuba

# Default goal of make file, this is built if you run just `make`
.DEFAULT_GOAL := all

# -- Define linker settings --
LDFLAGS := -L$(LHAPDF)/lib -L$(CUBA)/lib

# -- Define compiler settings --
# Setup warnings and optimizations
CXXFLAGS := -Wall -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter -O3 -std=c++11
# Include folders for libraries
CXXFLAGS := $(CXXFLAGS) -I$(LHAPDF)/include -I$(LOOPTOOLS)/include
CXXFLAGS := $(CXXFLAGS) -I$(CUBA)/include
# Include folders for project
CXXFLAGS := $(CXXFLAGS) -I./include

# -- Object file lists for the various target binaries -- 
# Splitting this off allows easier writing of utility tasks
#  and better dealing with dependency files (see common rules)

#general objects
#parameters
RESUM_OBJS := $(RESUM_OBJS) src/parameters.o src/Parameters_sm.o
#integration
RESUM_OBJS := $(RESUM_OBJS) src/MC/monte_carlo.o src/MC/4top_vegas.o
#PDFs
RESUM_OBJS := $(RESUM_OBJS) src/pdf/fit_coefficients.o src/pdf/deriv_pdf.o src/pdf/mellin_pdf.o 
#ulilities
RESUM_OBJS := $(RESUM_OBJS) src/utilities/inout.o src/utilities/read_slha.o
#hard functions
RESUM_OBJS := $(RESUM_OBJS) src/qqbar/qq_helamps.o src/qqbar/qq_process.o
#resumation
RESUM_OBJS := $(RESUM_OBJS) src/resum/resum_4top.o 

4top_program := $(RESUM_OBJS) 
4top_program := $(4top_program) programs/4top_run.o 

#PDFfits
PDF_program := $(RESUM_OBJS) programs/make_mellin_pdf.o


4top_program: $(4top_program)
	g++ -o ttH $(4top_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -lgfortran -lboost_program_options
	    	    
PDF: $(PDF_program)
	g++ -o PDF $(PDF_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -lgfortran -lboost_program_options
# -- Utility rules --

# all: builds all binaries
all: 4top_program

# clean: remove all generated files
clean:
	rm -f resum $(RESUM_OBJS) $(RESUM_OBJS:.o=.d)
	
# Mark the utility rules as phony to instruct make that they don't produce
#  files as output. This stops trouble when you accidentally create files
#  with those names.
.PHONY: all clean

# -- Rules for producing object files from source --
# Produce object files from c++ (the -MD ensures that this also creates
#  a dependency file for the source being compiled)
# The suffixes statement tells make that this is a general rule, and that
# .cpp and .o are file extensions.
.SUFFIXES: .cpp .o .f
.cpp.o:
	g++ -c $< -o $@ $(CXXFLAGS)

.f.o:
	gfortran -o $@ -c $<
# Include dependency files here to make changes to header files
#  cause proper recompilation of source files using them
#-include $(TRYLOOP_OBJS:.o=.d)
#-include $(RESUM_OBJS:.o=.d)
#-include $(FIXED_OBJS:.o=.d)
#-include $(MELLIN_OBJS:.o=.d)