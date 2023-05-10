ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTVERSION  := $(shell $(ROOTCONFIG) --version)
ROOTMAJORVERSION := $(word 1,$(subst ., ,$(ROOTVERSION)))
ROOTCINT=$(ROOTSYS)/bin/rootcint

# libraries generated in the current project
LIBDIR=$(ISR_UNFOLD_WD)/lib/

# header from V17 TUnfold package
HTUNFOLDV17=$(ISR_UNFOLD_WD)/TUnfold/

# sources for this project
SRC=$(ISR_UNFOLD_WD)/rootScripts/

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc -I$(HTUNFOLDV17) -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -fPIC $(ROOTCFLAGS)

LDFLAGS=$(ROOTLDFLAGS) -L$(LIBDIR) -Wl,-rpath $(LIBDIR)

LIB=unfold

_lib=$(LIBDIR)libisrunfold.so
lib:$(_lib)

SRCCODE=$(shell ls $(SRC)*.cpp)

_OBJ=$(SRCCODE:%.cpp=%_C.o)
OBJ=$(subst rootScripts,lib,$(_OBJ))

_DIC = $(SRCCODE:%.cpp=%_C_ACLiC_dict.cxx)                             
DIC = $(subst rootScripts,lib,$(_DIC))

# 
# make object files
#

$(LIBDIR)%_C_ACLiC_dict.cxx: $(SRC)%.h  $(SRC)Linkdef.h
	rootcling -f $@ -c -I$(HTUNFOLDV17) -p $^

$(LIBDIR)%_C.o: $(SRC)%.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(_lib): $(OBJ) $(DIC) 
	c++ $(CXXFLAGS) -shared -o $@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

clean:
	rm $(LIBDIR)ISR*
	rm $(LIBDIR)libisrunfold.so
