CXX = g++
CXXFLAGS = -g -Wall -std=c++11 -O3
vpath %.hpp ../include
DEPS = solver.hpp integrators.hpp system.hpp typedefs.hpp
INC = /home/fruzsina/oscode/include
NAGINC = /home/fruzsina/NAG/cll6i261dl/include
NAGSHARE = /home/fruzsina/NAG/cll6i261dl/lib/libnagc_nag.so

test-ms: test-ms.o test-main.o
	$(CXX) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(NAGSHARE)

test-airy: test-airy.o test-main.o
	$(CXX) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(NAGSHARE)

test-main.o: test-main.cpp
	$(CXX) -I$(INC) $(CXXFLAGS) -c -o $@ $<

test-airy.o: test-airy.cpp test-airy.hpp $(DEPS)
	$(CXX) -I$(NAGINC) $(CXXFLAGS) -c -o $@ $< 

test-ms.o: test-ms.cpp test-ms.hpp $(DEPS)
	$(CXX) -I$(NAGINC) -I$(INC) $(CXXFLAGS) -c -o $@ $< 

%.o: %.cpp %.hpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) *.o *~
