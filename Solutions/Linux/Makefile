SHELL = /bin/sh

export MAKE = make

.PHONY: default
default:
	@echo Please choose your platform from:
	@echo Type "make (win nix).(' '  clean install sos)".

#make all (add other platforms "Nix Nix Nix_double Win" ...)
PLAT := Nix
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

