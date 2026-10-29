SOURCES = bootstrapping.cpp contextchain.cpp seal.cpp extended_evaluator.cpp sparse_key_generator.cpp
SEAL_DIR = /home/paiclab/Documents/ckw/cic_hed/edit_distance_runner/seal
CXXFLAGS = -std=c++17 -w -I$(SEAL_DIR)/include/SEAL-4.1
HEADERS = $(wildcard *.h)

# The generated trace.cpp is one huge function, which gcc cannot compile for long traces.
# It is split into many small functions in SPLIT_DIR, which are compiled in parallel.
SPLIT_DIR = trace_split
LINES_PER_PART = 4000
JOBS ?= $(shell nproc)

bootstrapping_seal: $(SOURCES:.cpp=.o) constants.o $(SPLIT_DIR)/.stamp
	$(MAKE) -j$(JOBS) trace_objs
	g++ $(SOURCES:.cpp=.o) constants.o $(SPLIT_DIR)/*.o -L$(SEAL_DIR)/lib -lseal-4.1 -o bootstrapping_seal

$(SPLIT_DIR)/.stamp: trace.cpp split_trace.py
	rm -rf $(SPLIT_DIR)
	python3 split_trace.py trace.cpp $(SPLIT_DIR) $(LINES_PER_PART)
	touch $@

# Evaluated in the sub-make above, after the split files exist
TRACE_OBJS = $(patsubst %.cpp,%.o,$(wildcard $(SPLIT_DIR)/*.cpp))
trace_objs: $(TRACE_OBJS)

# The trace only calls into the bootstrapper, so heavy optimization has no measurable effect on timings
$(SPLIT_DIR)/%.o: $(SPLIT_DIR)/%.cpp $(HEADERS)
	ulimit -s unlimited; g++ $(CXXFLAGS) -I. -O1 -c $< -o $@

%.o: %.cpp $(HEADERS)
	g++ $(CXXFLAGS) -O3 -c $< -o $@

.PHONY: trace_objs clean distclean

# Keeps the library object files so that only the new trace needs to be compiled
clean:
	rm -f bootstrapping_seal
	rm -f constants.cpp constants.o
	rm -f tmp_trace_1.cpp
	rm -f tmp_trace_2.cpp
	rm -f tmp_trace_3.cpp
	rm -f tmp_trace_4.cpp
	rm -f trace.cpp trace.o
	rm -rf $(SPLIT_DIR)
	rm -f ops
	rm -f ctxt
	rm -f is_ntt_optimal_domain
	rm -r -f Files
	mkdir Files
	touch Files/.gitignore

distclean: clean
	rm -f *.o
