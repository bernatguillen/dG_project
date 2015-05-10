testobjects = unittest.o mesh1d.o
test = unittest


$(testobjects): CXXFLAGS = -g -Wall -L/usr/lib -lm -lgsl
$(test): CXXFLAGS = -g -Wall -L/usr/lib -lm


.PHONY: all clean depend

all: $(test)
	echo make complete

$(test) : $(testobjects)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgsl


clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) $(jactargets) $(multargets)
	$(RM) *.dat
depend:

	$(CXX) -MM $(CXXFLAGS) *.cc > .depend


-include .depend
