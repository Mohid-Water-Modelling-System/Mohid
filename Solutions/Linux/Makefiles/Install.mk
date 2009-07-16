SOURCE := $(TARGET)

.PHONY: install

install: $(SOURCE)
	@-$(CP) $< $(DESTDIR)/`date +%G%m%d`_$(addprefix $(VER), $<)
	@echo Installed $<.
