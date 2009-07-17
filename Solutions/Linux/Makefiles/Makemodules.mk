#make all
MODALL = $(addsuffix .all, $(MODULES))
.PHONY: $(MODALL)
$(MODALL):
	@-$(MAKE) -C $(@:.all=) all

.PHONY: all
all: $(MODALL)

#make clean
MODCLEAN = $(addsuffix .clean, $(MODULES))
.PHONY: $(MODCLEAN)
$(MODCLEAN):
	@-$(MAKE) -C $(@:.clean=) clean

.PHONY: clean
clean: $(MODCLEAN)

#make install
MODINSTALL = $(addsuffix .install, $(MODULES))
.PHONY: $(MODINSTALL)
$(MODINSTALL) :
	@-$(MAKE) -C $(@:.install=) install

.PHONY: install
install: $(MODINSTALL)
