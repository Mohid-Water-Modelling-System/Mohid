OBJS = $(FILES:.$(S)=.$(O))
SRCF = $(SRCDIR)

.PHONY: all

all: $(TARGET)

$(TARGET) : $(OBJS) $(LIBS)
	$(CC) $(LFLAGS) -o $@ $^ $(LLFLAGS)
	@echo Finished building $@.

%.$(O) : $(SRCF)/%.$(F)
	$(CC) $(CCFLAGS) $(INCS) $<
	@echo $* .................. [OK]

