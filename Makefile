testobjects = unittest.o mesh1d.o
test = unittest
advexampleobjects = mesh1d.o d1Advection.o Advection_example.o
advexample = Adv_1D
objects = $(testobjects) $(advexampleobjects)
targets = $(test) $(advexample)

$(objects): CXXFLAGS = -g -Wall -L/usr/lib -lm -lgsl
$(targets): CXXFLAGS = -g -Wall -L/usr/lib -lm


.PHONY: all clean depend

all: $(targets)
	echo make complete

$(test) : $(testobjects)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgsl

$(advexample) : $(advexampleobjects)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgsl

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) $(jactargets) $(multargets)
	$(RM) *.dat
	$(RM) *.stackdump
depend:

	$(CXX) -MM $(CXXFLAGS) *.cc > .depend


-include .depend
