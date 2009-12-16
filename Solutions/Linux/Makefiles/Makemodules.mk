#make all
MODALL = $(addsuffix .all, $(MODULES))
.PHONY: $(MODALL)
$(MODALL):
	@-$(MAKE) $(MKFLAGS) -C $(@:.all=) all

.PHONY: all
all: $(MODALL)

#make clean
MODCLEAN = $(addsuffix .clean, $(MODULES))
.PHONY: $(MODCLEAN)
$(MODCLEAN):
	@-$(MAKE) $(MKFLAGS) -C $(@:.clean=) clean

.PHONY: clean
clean: $(MODCLEAN)

#make install
MODINSTALL = $(addsuffix .install, $(MODULES))
.PHONY: $(MODINSTALL)
$(MODINSTALL) :
	@-$(MAKE) $(MKFLAGS) -C $(@:.install=) install

.PHONY: install
install: $(MODINSTALL)
