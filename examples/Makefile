CXX = g++
CXXFLAGS = -g -std=c++11 -O3 #-Wall
inc = include
deps = ../include/system.hpp ../include/rksolver.hpp ../include/wkbsolver.hpp ../include/interpolator.hpp ../include/solver.hpp

vdp: vdp.o
	$(CXX) -I../include/ $(CXXFLAGS) -o $@ $^

burst: burst.o
	$(CXX) -I../include/ $(CXXFLAGS) -o $@ $^

issue13: issue13.o
	$(CXX) -I../include/ $(CXXFLAGS) -o $@ $^

issue15: issue15.o
	$(CXX) -I../include/ $(CXXFLAGS) -o $@ $^

airy: airy.o
	$(CXX) -I../include/ $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.hpp $(deps)
	$(CXX) -I../include/ $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp $(deps)
	$(CXX) -I../include/ $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) *.o
