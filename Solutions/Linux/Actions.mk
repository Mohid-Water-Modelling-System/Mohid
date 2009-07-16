#make all
MODALL = $(addsuffix .all, $(MODULES))
.PHONY: all $(MODALL)
all: $(MODALL)
$(MODALL):
	-$(MAKE) -C $(@:.all=) all

#make clean
MODCLEAN = $(addsuffix .clean, $(MODULES))
.PHONY: clean $(MODCLEAN)
clean: $(MODCLEAN)
$(MODCLEAN):
	-$(MAKE) -C $(@:.clean=) clean

#make install
MODINSTALL = $(addsuffix .install, $(MODULES))
.PHONY: install $(MODINSTALL)
install: $(MODINSTALL)
$(MODINSTALL) :
	-$(MAKE) -C $(@:.install=) install

