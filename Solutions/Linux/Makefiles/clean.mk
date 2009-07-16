.PHONY: clean
clean:
	@-$(DEL) *.$(O) *.$(MOD) $(TARGET)
	@echo erased $(TARGET) files.

