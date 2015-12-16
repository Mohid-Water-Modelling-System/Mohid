SHELL = /bin/sh

export MAKE = make

.PHONY: default
default:
	@echo Please choose your platform from:
	@echo Type "make Nix.(' '  clean install)".

#make all 
PLAT := Nix
.PHONY : $(PLAT)
$(PLAT):
	@$(MAKE) -f $@.mk all

#make clean
PLATCLEAN := $(addsuffix .clean, $(PLAT))
.PHONY : $(PLATCLEAN)
$(PLATCLEAN) :
	@$(MAKE) -f $(@:.clean=).mk clean

#make install
PLATINST := $(addsuffix .install, $(PLAT))
.PHONY : $(PLATINST)
$(PLATINST) :
	@$(MAKE) -f $(@:.install=).mk install

