CXX = g++
CXXFLAGS = -g -Wall -std=c++11 -O3
#inc = include
deps =  system.hpp rksolver.hpp wkbsolver.hpp interpolator.hpp

solver: solver.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.hpp $(deps)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) $(test)/*.o
