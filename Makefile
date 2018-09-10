CXX = g++
CXXFLAGS = -g -Wall -std=c++11 -O3
VPATH = include/
DEPS = solver.hpp integrators.hpp system.hpp typedefs.hpp
INCDIR = /home/fruzsina/NAG/cll6i261dl/include
SHARDIR = /home/fruzsina/NAG/cll6i261dl/lib/libnagc_nag.so

nag_rk_minimal: nag_rk_minimal.cpp
	$(CXX) -I$(INCDIR) $< $(SHARDIR) -std=c++11 -lpthread -lm -o $@

test_nag: test_nag.cpp
	$(CXX) -I$(INCDIR) $< $(SHARDIR) -std=c++11 -lpthread -lm -o $@

test-airy: test-airy.o test-main.o
	$(CXX) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(SHARDIR)

test-main.o: test-main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test-airy.o: test-airy.cpp test-airy.hpp $(DEPS)
	$(CXX) -I$(INCDIR) $(CXXFLAGS) -c -o $@ $< 

%.o: %.cpp %.hpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o *~
