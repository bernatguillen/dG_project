testobjects = unittest.o mesh1d.o
test = unittest


$(objects): CXXFLAGS = -g -Wall -lm
$(targets): CXXFLAGS = -g -Wall -lm


.PHONY: all clean depend

all: $(test)
	echo make complete

$(test) : $(testobjects)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) $(jactargets) $(multargets)
	$(RM) *.dat
depend:

	$(CXX) -MM $(CXXFLAGS) *.cc > .depend


-include .depend
