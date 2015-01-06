OBJS1 = $(FILES:.$(S)=.$(Obj))
OBJS = $(OBJS1:.F=.$(Obj)) 
SRCF = $(SRCDIR)
SUFF := $(suffix $(TARGET))

#Make all
.PHONY: all
all: $(TARGET)

$(TARGET) : $(OBJS) $(LIBS)
	@echo  
	@echo Build with compilation flags $(CCFLAGS) $(INCS)
	@echo 
ifeq ($(SUFF),$(SUFFLIB))
	$(AR) $@ $^
else
	$(CC) $(LFLAGS) -o $@ $^ $(LLFLAGS)
endif
	@echo 
	@echo Finished building $@.
	@echo 

#Fortran compilation rule
%.$(Obj) : $(SRCF)/%.$(F)
	@$(CC) $(CCFLAGS) $(INCS) $<
	@echo $* .................. [OK]

#make clean
.PHONY: clean
clean:
	@-$(DEL) *.$(Obj) *.$(MOD) $(TARGET)
	@echo erased $(TARGET) files.

#make install
SOURCE := $(TARGET)
.PHONY: install
install: $(SOURCE)
	@-$(CP) $< $(DESTDIR)/`date +%G%m%d`_$(addprefix $(VER), $<)
	@echo Installed $<.
