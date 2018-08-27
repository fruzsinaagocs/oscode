CXX = g++ 
CXXFLAGS = -g -Wall
VPATH = include/
DEPS = solver.hpp integrators.hpp system.hpp typedefs.hpp

test-airy: test-airy.o test-main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test-main.o: test-main.cpp
	$(CXX)$ -c -o $@ $<

#test: test.o
#	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.hpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o *~
