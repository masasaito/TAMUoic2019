
INC_COMMON = ./Common.mk
include $(INC_COMMON)

## Macros
TARGETS  = tamuoic2019encode
SHLS     = $(STARGETS)
BUPIDX   = $(shell date "+%Y%m%d_%H%M")
BUPFILE  = tamuoic2019-src.$(BUPIDX).tar.gz

## Files
SRCS = main_TAMUoic2019_encode.f90 
OBJS = $(SRCS:.f90=.o)

## Commands
all : $(TARGETS) $(STARGETS)

install : $(TARGETS) $(STARGETS)
	@ mkdir -p $(BINDIR)
	cp $(TARGETS) $(STARGETS) $(BINDIR)
clean  :
	rm -f *.o *.mod $(TARGETS) *~
echo   :
	@ echo $(SRCS) $(STARGETS)
print  :
	@ cat $(SRCS) $(STARGETS)
backup : $(BUPFILE)
$(BUPFILE) :
	@ rm -f $(BUPFILE)
	tar cvzf $(BUPFILE) $(SRCS) $(STARGETS) Makefile *.mk

## Suffix rules
.f90.o :
	$(FC) -c $(FCFLAGS) $(LIB_FLAGS) $<

## Targets
tamuoic2019encode    : main_TAMUoic2019_encode.o  $(LIB_OBJS) 
	$(FC) -o $@ $(FCFLAGS) $^

## Dependency
main_TAMUoic2019_encode.o         : $(LIB_OBJS)
