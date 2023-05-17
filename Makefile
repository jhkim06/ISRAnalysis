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

# headers from V17 TUnfold package
HTUNFOLDV17=$(ISR_UNFOLD_WD)/TUnfold/

# sources for this project
SRC=$(ISR_UNFOLD_WD)/src/
INCLUDE=$(ISR_UNFOLD_WD)/include/

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(HTUNFOLDV17) -I$(INCLUDE) -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -fPIC $(ROOTCFLAGS)
LDFLAGS=$(ROOTLDFLAGS) -L$(LIBDIR) -Wl,-rpath $(LIBDIR)
LIB=unfold

SRCCODE=$(shell ls $(SRC)*.cpp)

_lib=$(LIBDIR)libisrunfold.so
lib:$(_lib)

_OBJ=$(SRCCODE:%.cpp=%_C.o)
OBJ=$(subst src,lib,$(_OBJ))

_DIC = $(SRCCODE:%.cpp=%_C_ACLiC_dict.cxx)                             
DIC = $(subst src,lib,$(_DIC))

# make object files

$(LIBDIR)%_C.o: $(SRC)%.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIBDIR)%_C_ACLiC_dict.cxx: $(INCLUDE)%.h  $(INCLUDE)Linkdef.h
	rootcling -f $@ -c -I$(HTUNFOLDV17) -p $^

$(_lib): $(OBJ) $(DIC) 
	c++ $(CXXFLAGS) -shared -o $@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

clean:
	rm $(LIBDIR)ISR*
	rm $(LIBDIR)libisrunfold.so
