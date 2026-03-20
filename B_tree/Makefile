NODES := nodes
BENCHMARK := benchmark
TEST := test
BUILD := build
OBJECTS := objects

ALL_LDFLAGS := $(LDFLAGS) $(EXTRA_LDFLAGS)

clean:
	@rm -rf $(OBJECTS) $(BUILD) $(NODES)

# Directories
$(NODES) $(BENCHMARK) $(TEST) $(BUILD) $(OBJECTS):
	$(Q)mkdir -p $@

$(OBJECTS)/btree.o: btree.cpp btree.h | $(OBJECTS) $(NODES)
	@$(CXX) $(CXXFLAGS) -g -fPIC -c btree.cpp -o $@

$(BUILD)/libbtree.so: $(OBJECTS)/btree.o | $(BUILD)
	@$(CXX) $(CXXFLAGS) -shared  -fPIC -g $^ -o $@	

$(OBJECTS)/main.o: main.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c main.cpp -o $@

$(BUILD)/main: $(OBJECTS)/main.o $(BUILD)/libbtree.so | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $(OBJECTS)/main.o $(ALL_LDFLAGS) -g -lelf -lz -L$(BUILD) -lbtree -Wl,-rpath,. -o $@ 

$(OBJECTS)/test.o: $(TEST)/test.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c $(TEST)/test.cpp -o $@

$(BUILD)/test: $(OBJECTS)/test.o $(BUILD)/libbtree.so | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $(OBJECTS)/test.o $(ALL_LDFLAGS) -g -L$(BUILD) -lbtree -Wl,-rpath,. -lgtest -lgmock -pthread -lelf -lz -o $@

$(OBJECTS)/benchmark.o: $(BENCHMARK)/benchmark.cpp btree.h | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -g -c $(BENCHMARK)/benchmark.cpp -o $@

$(BUILD)/benchmark: $(OBJECTS)/benchmark.o $(BUILD)/libbtree.so | $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $(OBJECTS)/benchmark.o $(ALL_LDFLAGS) -g -std=c++14  -L$(BUILD) -lbtree -Wl,-rpath,. -lbenchmark -lpthread -lelf -lz -o $@

complete: $(BUILD)/main $(BUILD)/test $(BUILD)/benchmark
	@rm -rf $(OBJECTS)
	cp btree.h $(BUILD)

.DELETE_ON_ERROR:

.SECONDARY:
