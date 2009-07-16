OBJS = $(FILES:.$(S)=.$(O))
SRCF = $(SRCDIR)

.PHONY: all

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(AR) $@ $^
	@echo Finished building $@.

%.$(O) : $(SRCF)/%.$(F)
	$(CC) $(CCFLAGS) $(INCS) $<
	@echo $* .................. [OK]

