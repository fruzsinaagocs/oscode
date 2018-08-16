CXX = g++ 
CXXFLAGS = -g -Wall
DEPS = solver.hpp integrators.hpp system.hpp typedefs.hpp

test: test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.hpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o *~
