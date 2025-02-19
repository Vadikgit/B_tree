OBJECTS := objects
NODES := nodes
ALL_LDFLAGS := $(LDFLAGS) $(EXTRA_LDFLAGS)

clean:
	@rm -rf $(OBJECTS)/* main test benchmark
	@rm -rf $(OBJECTS)
	@rm -rf $(NODES)/*
	@rm -rf $(NODES)

# Directories
$(OBJECTS) $(NODES):
	$(Q)mkdir -p $@

$(OBJECTS)/btree.o: btree.cpp btree.h | $(OBJECTS) $(NODES)
	@$(CXX) $(CXXFLAGS) -g -c btree.cpp -o $@

$(OBJECTS)/main.o: main.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c main.cpp -o $@

$(OBJECTS)/test.o: test.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c test.cpp -o $@
	
$(OBJECTS)/benchmark.o: benchmark.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c benchmark.cpp -o $@

main: $(OBJECTS)/btree.o $(OBJECTS)/main.o | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $^ $(ALL_LDFLAGS) -g -lelf -lz -o $@

test: $(OBJECTS)/btree.o $(OBJECTS)/test.o | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $^ $(ALL_LDFLAGS) -g  -lgtest -lgmock -pthread -lelf -lz -o $@

benchmark: $(OBJECTS)/btree.o $(OBJECTS)/benchmark.o | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $^ $(ALL_LDFLAGS) -g -std=c++14 -lbenchmark -lpthread -lelf -lz -o $@

.DELETE_ON_ERROR:

.SECONDARY:
