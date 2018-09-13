CXX = g++
CXXFLAGS = -g -Wall -std=c++11 -O3
test = test
inc = include
naginc = /home/fruzsina/NAG/cll6i261dl/include
nagshare = /home/fruzsina/NAG/cll6i261dl/lib/libnagc_nag.so
deps = $(inc)/solver.hpp $(inc)/integrators.hpp $(inc)/system.hpp $(inc)/typedefs.hpp

$(test)/%: $(test)/%.o $(test)/test-main.o
	$(CXX) -I$(inc) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(nagshare)

$(test)/time-ms: $(test)/time-ms.o
	$(CXX) -I$(inc) $(CXXFLAGS) -lpthread -lm -o $@ $^ $(nagshare)

$(test)/test-main.o: $(test)/test-main.cpp
	$(CXX) -I$(inc) $(CXXFLAGS) -c -o $@ $<

$(test)/time-ms.o: $(test)/time-ms.cpp $(test)/test-ms.hpp $(deps)
	$(CXX) -I$(naginc) -I$(inc) $(CXXFLAGS) -c -o $@ $<

$(test)/%.o: $(test)/%.cpp $(test)/%.hpp $(deps)
	$(CXX) -I$(naginc) -I$(inc) $(CXXFLAGS) -c -o $@ $<

clean: 
	$(RM) $(test)/*.o $(test)/*~
