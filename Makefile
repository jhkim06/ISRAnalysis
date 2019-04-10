ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTVERSION  := $(shell $(ROOTCONFIG) --version)
ROOTMAJORVERSION := $(word 1,$(subst ., ,$(ROOTVERSION)))
ROOTCINT=$(ROOTSYS)/bin/rootcint

# libraries generated in the current project
LIBDIR=/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/

# header from V17 TUnfold package
HTUNFOLDV17=/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/TUnfold/

# sources for this project
SRC=/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc -I$(HTUNFOLDV17) -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -fPIC $(ROOTCFLAGS)

LDFLAGS=$(ROOTLDFLAGS) -L$(LIBDIR) -Wl,-rpath $(LIBDIR)

LIB=unfold

SRCCODE=$(shell ls $(SRC)ISR_*.C)

_OBJ=$(SRCCODE:%.C=%_C.o)
OBJ=$(subst rootScripts,lib,$(_OBJ))

_DIC = $(SRCCODE:%.C=%_C_ACLiC_dict.cxx)                             
DIC = $(subst rootScripts,lib,$(_DIC))

# 
# make object files
#

$(LIBDIR)%_C_ACLiC_dict.cxx: $(SRC)%.h $(SRC)Linkdef.h
	rootcling -f $@ -c `root-config --ldflags` -I$(HTUNFOLDV17) -p $^

$(LIBDIR)%_C.o: $(SRC)%.C 
	$(CXX) $(CXXFLAGS) -c $< -o $@

libisrunfold.so: $(OBJ) $(DIC) 
	c++ $(CXXFLAGS) -shared -o $(LIBDIR)$@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

clean:
	rm $(LIBDIR)ISR_*
	rm $(LIBDIR)libisrunfold.so
