
#========================#
# check env has been set #
#========================#
ifeq ($(strip $(DJANGO_CPP)),)
$(info It seems DjanGo environment has not been set.)
$(info Source one of the env files in ../../../env or create a new one for your machine.)
$(error DJANGO_CPP is set to $(DJANGO_CPP))
#else
#$(info Django will be compiled with $(DJANGO_CPP))
endif

#=====================#
#     directories     #
#=====================#

# tool sources
SRC_DIR        = ../src

# EX2D sources 
EX2D_DIR       = ../../../misc/external/EX2D

# binaries
BIN_DIR        = ../../../bin

# include
DJANGO_INC     = ../../../include
NR3_INC        = ../../../misc/external/Numerical-Recipes

#=====================#
#       vpaths        #
#=====================#

VPATH = $(SRC_DIR):$(EX2D_DIR)

#=====================#
#       options       #
#=====================#

# SINGLE / DOUBLE PRECISION
ifeq ($(findstring double,$(precision)), double)
DJANGO_CPP_FLAGCOMP+=-D_DOUBLE_PRECISION_
DJANGO_CUDA_FLAGCOMP+=-D_DOUBLE_PRECISION_
DJANGO_HIP_FLAGCOMP+=-D_DOUBLE_PRECISION_
BUILD_MESSAGE='*** BUILD OK (DOUBLE PRECISION) *** '
else
BUILD_MESSAGE='*** BUILD OK (SINGLE PRECISION) *** '
endif

# WITH CUDA
ifneq ($(strip $(DJANGO_CUDA)),)
DJANGO_CPP_FLAGCOMP+=-D__CUDA__
endif

# WITH DPC++
ifeq ($(findstring dpcpp,$(DJANGO_CPP)), dpcpp)
DJANGO_CPP_FLAGCOMP+=-D__DPCPP__
endif

# WITH HIP
ifneq ($(strip $(DJANGO_HIP)),)
DJANGO_CPP_FLAGCOMP+=-D__HIP__
endif

# WITH OPENACC
ifneq ($(strip $(DJANGO_CPP_OPENACC_FLAG)),)
DJANGO_CPP_FLAGCOMP+=-D__OPENACC__
DJANGO_CPP_FLAGCOMP+=$(DJANGO_CPP_OPENACC_FLAG)
endif

#$(info DJANGO_CPP_FLAGCOMP="$(DJANGO_CPP_FLAGCOMP)")

DJANGO_CPP_FLAGCOMP+= -I ${DJANGO_INC} -I$(NR3_INC) $(COPTION)

#=====================#
#     objects         #
#=====================#

#.SUFFIXES: .o .cpp
.cpp.o:
	$(DJANGO_CPP) -c $(DJANGO_CPP_FLAGCOMP) $<
.f.o:
	$(DJANGO_FC) -c $(DJANGO_CPP_FLAGCOMP) $<

all: analytic_homo analytic_layer build_acqui_file build_acqui_file_surround \
build_model_gradient build_model_homo build_model_sphere build_model_layer \
check_file check_rms compute_element_matrices \
convert_1D_to_2D double2single eigen_sol \
mulgrid single2double subgrid transpose transpose_2d \
wavelet EX2DDIR_Berg_If EX2DVAEL_Berg_If EX2DELEL_Berg_If
	@echo $(BUILD_MESSAGE)

EX2DDIR_Berg_If: EX2DDIR_Berg_If.o
	@mkdir -p $(BIN_DIR)
	$(DJANGO_FC) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ EX2DDIR_Berg_If.o $(LIBS)

EX2DVAEL_Berg_If: EX2DVAEL_Berg_If.o
	$(DJANGO_FC) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ EX2DVAEL_Berg_If.o $(LIBS)

EX2DELEL_Berg_If: EX2DELEL_Berg_If.o
	$(DJANGO_FC) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ EX2DELEL_Berg_If.o $(LIBS)

analytic_homo: analytic_homo.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ analytic_homo.o $(LIBS)

analytic_layer: analytic_layer.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ analytic_layer.o $(LIBS)

build_acqui_file: build_acqui_file.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_acqui_file.o $(LIBS)

build_acqui_file_surround: build_acqui_file_surround.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_acqui_file_surround.o $(LIBS)

build_model_gradient: build_model_gradient.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_model_gradient.o $(LIBS)

build_model_homo: build_model_homo.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_model_homo.o $(LIBS)

build_model_layer: build_model_layer.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_model_layer.o $(LIBS)

build_model_sphere: build_model_sphere.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ build_model_sphere.o $(LIBS)

check_file: check_file.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ check_file.o $(LIBS)

check_rms: check_rms.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ check_rms.o $(LIBS)

compute_element_matrices: compute_element_matrices.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ compute_element_matrices.o $(LIBS)

convert_1D_to_2D: convert_1D_to_2D.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ convert_1D_to_2D.o $(LIBS)

double2single: double2single.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ double2single.o $(LIBS)

eigen_sol: eigen_sol.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ eigen_sol.o $(LIBS)

mulgrid: mulgrid.o
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ mulgrid.o $(LIBS)

single2double: single2double.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ single2double.o $(LIBS)

subgrid: subgrid.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ subgrid.o $(LIBS)

transpose: transpose.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ transpose.o $(LIBS)

transpose_2d: transpose_2d.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ transpose_2d.o $(LIBS)

wavelet: wavelet.o 
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ wavelet.o $(LIBS)

clean:
	rm -f *.o *.mod *~ 







