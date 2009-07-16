SHELL = /bin/sh

export MAKE = make

.PHONY: default
default:
	@echo Please choose your platform from:
	@echo Type "make (win nix).(' '  clean install sos)".

#make all
PLAT := nix win nix_double nix_mpi nix_openmp nix_pgf nix_fast nix_single
.PHONY : $(PLAT)
$(PLAT):
	$(MAKE) -f $@.mk all

#make clean
PLATCLEAN := $(addsuffix .clean, $(PLAT))
.PHONY : $(PLATCLEAN)
$(PLATCLEAN) :
	$(MAKE) -f $(@:.clean=).mk clean

#make install
PLATINST := $(addsuffix .install, $(PLAT))
.PHONY : $(PLATINST)
$(PLATINST) :
	$(MAKE) -f $(@:.install=).mk install

