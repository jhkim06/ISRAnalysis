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

# use % later to shorten the length
# make object files
ISR_histTUnfold_C.o: $(SRC)ISR_histTUnfold.C
	c++  $(CXXFLAGS) -c $^ -o $(LIBDIR)$@

ISR_saveHists_C.o: $(SRC)ISR_saveHists.C
	c++  $(CXXFLAGS) -c $^ -o $(LIBDIR)$@

ISR_unfoldUtils_C.o: $(SRC)ISR_unfoldUtils.C
	c++  $(CXXFLAGS) -c $^ -o $(LIBDIR)$@

ISR_drawUtils_C.o: $(SRC)ISR_drawUtils.C
	c++  $(CXXFLAGS) -c $^ -o $(LIBDIR)$@

# make dictionary
ISR_histTUnfold_C_ACLiC_dict.cxx: $(SRC)ISR_histTUnfold.h $(SRC)Linkdef.h
	rootcling -f $(LIBDIR)$@ -c `root-config --ldflags`  -I$(HTUNFOLDV17) -p $^

ISR_saveHists_C_ACLiC_dict.cxx: $(SRC)ISR_saveHists.h $(SRC)Linkdef.h
	rootcling -f $(LIBDIR)$@ -c `root-config --ldflags`  -I$(HTUNFOLDV17) -p $^

ISR_unfoldUtils_C_ACLiC_dict.cxx: $(SRC)ISR_unfoldUtils.h $(SRC)Linkdef.h
	rootcling -f $(LIBDIR)$@ -c `root-config --ldflags`  -I$(HTUNFOLDV17) -p $^

ISR_drawUtils_C_ACLiC_dict.cxx: $(SRC)ISR_drawUtils.h $(SRC)Linkdef.h
	rootcling -f $(LIBDIR)$@ -c `root-config --ldflags`  -I$(HTUNFOLDV17) -p $^

# make shared library to be loaded in python 
libhistTUnfold_C.so: $(LIBDIR)ISR_histTUnfold_C_ACLiC_dict.cxx $(LIBDIR)ISR_histTUnfold_C.o
	c++ $(CXXFLAGS) -shared -o $(LIBDIR)$@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

libsaveHists_C.so: $(LIBDIR)ISR_saveHists_C_ACLiC_dict.cxx $(LIBDIR)ISR_saveHists_C.o
	c++ $(CXXFLAGS) -shared -o $(LIBDIR)$@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

libunfoldUtils_C.so: $(LIBDIR)ISR_unfoldUtils_C_ACLiC_dict.cxx $(LIBDIR)ISR_unfoldUtils_C.o
	c++ $(CXXFLAGS) -shared -o $(LIBDIR)$@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)

libdrawUtils_C.so: $(LIBDIR)ISR_drawUtils_C_ACLiC_dict.cxx $(LIBDIR)ISR_drawUtils_C.o
	c++ $(CXXFLAGS) -shared -o $(LIBDIR)$@ $^ $(LDFLAGS) -l$(LIB) \
	$(ROOTLIBS)


clean:
	rm $(LIBDIR)/ISR_*
	rm $(LIBDIR)/libhistTUnfold*
