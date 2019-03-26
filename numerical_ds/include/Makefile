CXX = mpiCC
CXXFLAGS = -g -Wall -std=c++11 -O3
naginc = /home/fruzsina/NAG/cll6i261dl/include
nagshare = /home/fruzsina/NAG/cll6i261dl/lib/libnagc_nag.so
deps =  system.hpp rksolver.hpp wkbsolver.hpp interpolator.hpp solver.hpp

bingo-%: bingo-%.o
	$(CXX) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(nagshare)

ms-solver-%: ms-solver-%.o
	$(CXX) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(nagshare)

solver: solver.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.hpp $(deps)
	$(CXX) -I$(naginc) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp $(deps)
	$(CXX) -I$(naginc) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o
