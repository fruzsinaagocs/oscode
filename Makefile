CXX = g++
CXXFLAGS = -g -Wall -std=c++11
VPATH = include/
DEPS = solver.hpp integrators.hpp system.hpp typedefs.hpp

test-airy: test-airy.o test-main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test-main.o: test-main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp %.hpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o *~
