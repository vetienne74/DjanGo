
#========================#
# check env has been set #
#========================#
ifeq ($(strip $(DJANGO_CPP)),)
$(info It seems DjanGo environment has not been set.)
$(info Source one of the env files in ../env or create a new one for your machine.)
$(error DJANGO_CPP is set to $(DJANGO_CPP))
#else
#$(info Django will be compiled with $(DJANGO_CPP))
endif

#=====================#
#     directories     #
#=====================#

# source files
SRC_DIR        = ../src
EXPAT_SRC_DIR  = ../misc/external/expat-2-2-0

# header files
INC_DIR        = ../include
EXPAT_INC_DIR  = ../misc/external/expat-2-2-0
NR3_INC_DIR    = ../misc/external/Numerical-Recipes

# executable files
BIN_DIR        = ../bin

#=====================#
#       vpaths        #
#=====================#

VPATH = $(SRC_DIR):$(EXPAT_SRC_DIR)

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

DJANGO_CPP_FLAGCOMP+=-I. -I$(INC_DIR) -I$(EXPAT_INC_DIR) -I$(NR3_INC_DIR) -I${SRC_DIR} $(COPTION) -DHAVE_EXPAT_CONFIG_H

DJANGO_GIT_VER=1

#=====================#
#     objects         #
#=====================#

# all objects
CPP_OBJ = program.o program_modelling.o program_fwi.o program_rtm.o program_guitar.o \
scheme.o scheme_factory.o domain.o \
singleton.o output_report.o inversion.o main.o \
parse_argument.o \
modelling.o \
FDM.o \
FDM_1D.o \
FDM_1D_ac_iso.o \
FDM_1D_ac_iso_1st.o FDM_1D_ac_iso_1st_schemes.o \
FDM_1D_ac_iso_2nd.o FDM_1D_ac_iso_2nd_schemes.o \
FDM_1D_ac_lossy.o \
FDM_1D_ac_lossy_1st.o FDM_1D_ac_lossy_1st_schemes.o \
FDM_1D_ac_lossy_2nd.o FDM_1D_ac_lossy_2nd_schemes.o \
FDM_2D.o \
FDM_2D_ac_iso.o \
FDM_2D_ac_iso_1st.o FDM_2D_ac_iso_1st_schemes.o \
FDM_2D_ac_iso_2nd.o FDM_2D_ac_iso_2nd_schemes.o \
FDM_2D_ac_lossy.o \
FDM_2D_ac_lossy_1st.o FDM_2D_ac_lossy_1st_schemes.o \
FDM_2D_ac_lossy_2nd.o FDM_2D_ac_lossy_2nd_schemes.o \
FDM_2D_el_iso.o \
FDM_2D_el_iso_1st.o FDM_2D_el_iso_1st_schemes.o \
FEM.o \
FEM_1D.o \
FEM_1D_1st.o FEM_1D_1st_ac_iso.o FEM_1D_1st_ac_lossy.o \
FEM_1D_2nd.o FEM_1D_2nd_ac_iso.o \
FEM_2D.o \
FEM_2D_1st.o FEM_2D_1st_ac_iso.o FEM_2D_1st_ac_lossy.o \
FEM_2D_1st_el_iso.o \
FEM_2D_2nd.o FEM_2D_2nd_ac_iso.o \
acquisition.o \
freq_group.o \
grid.o \
grid_1D.o grid_1D_int.o grid_1D_float.o grid_1D_complex.o \
grid_2D.o grid_2D_int.o grid_2D_float.o grid_2D_complex.o \
grid_3D.o grid_3D_float.o grid_3D_complex.o \
model.o \
data.o data_std.o \
gradient.o \
parse_xml.o xmlparse.o xmlrole.o xmltok.o \
boundary.o variable.o snapshot.o

# build list of needed objects
#-----------------------------
# common C++ objects
ALL_OBJ = ${CPP_OBJ}

.SUFFIXES: .o .cpp .c

.cpp.o:
	$(DJANGO_CPP) -c $(DJANGO_CPP_FLAGCOMP) $<
.c.o:
	$(DJANGO_CPP) -c $(DJANGO_CPP_FLAGCOMP) $<

all: version django

django: $(ALL_OBJ)
	@mkdir -p $(BIN_DIR)
	$(DJANGO_CPP) $(DJANGO_CPP_FLAGCOMP) -o $(BIN_DIR)/$@ $(ALL_OBJ) $(DJANGO_CPP_LIB)
	@echo $(BUILD_MESSAGE)

clean:
	rm -f *.o *.mod *~ ./version $(BIN_DIR)/django

version:
	@if [ ${DJANGO_GIT_VER} ] ; then \
		./git_info.sh ../ > $@ ;\
		if test ! -e version ; then \
			echo "No version found, initialization" ;\
			cp version $(SRC_DIR)/version_django.h ; \
		fi;\
		if ! diff -q version $(SRC_DIR)/version_django.h ; then \
			echo "version has changed, regenerating version info"; \
			cp version $(SRC_DIR)/version_django.h ; \
		else \
			echo "version info is up to date"; \
		fi;\
	fi







