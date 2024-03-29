################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = test_advect1d
OBJS = 

LIBS = $(LIBDIR)/libScaleFECore.a

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) makebin

makedir:
	mkdir -p $(BUILD_DIR)

makebin: $(BINNAME)
	@echo "Complete making."

run:
	export OMP_NUM_THREADS=1
	mpirun -n 1 ./$(BINNAME) advect.conf

vis:
	bash ./visualize/visualize.sh

make_sample_aux:
	make -C $(SAMPLE_AUX_DIR)

$(BINNAME): $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o: $(BINNAME).f90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

distclean: clean
	rm -f $(BINNAME)

clean:
	rm -rf $(BUILD_DIR)

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.mod: %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include -I$(SAMPLE_AUX_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments
