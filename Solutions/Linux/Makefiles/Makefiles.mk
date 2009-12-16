OBJS1 = $(FILES:.$(S)=.$(O))
OBJS = $(OBJS1:.F=.$(O)) 
SRCF = $(SRCDIR)
SUFF := $(suffix $(TARGET))

#Make all
.PHONY: all
all: $(TARGET)

$(TARGET) : $(OBJS) $(LIBS)
ifeq ($(SUFF),$(SUFFLIB))
	@$(AR) $@ $^
else
	@$(CC) $(LFLAGS) -o $@ $^ $(LLFLAGS)
endif
	@echo Finished building $@.

#Fortran compilation rule
%.$(O) : $(SRCF)/%.$(F)
	@$(CC) $(CCFLAGS) $(INCS) $<
	@echo $* .................. [OK]

#make clean
.PHONY: clean
clean:
	@-$(DEL) *.$(O) *.$(MOD) $(TARGET)
	@echo erased $(TARGET) files.

#make install
SOURCE := $(TARGET)
.PHONY: install
install: $(SOURCE)
	@-$(CP) $< $(DESTDIR)/`date +%G%m%d`_$(addprefix $(VER), $<)
	@echo Installed $<.
