include Makefile.arch

SOURCES=$(wildcard *.cc) $(wildcard MT2/*.cc)

OBJECTS=$(SOURCES:.cc=.o)

LIB=libTools.so

$(LIB): $(OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -o $@

%.o:  %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: $(LIB) 
clean:
	rm -f *.o \ 
	rm -f *.d \
	rm -f *.so \
