CXX = g++
CXXFLAGS = -O3

SRC = src/main.cpp src/reader.cpp src/compare.cpp src/tree.cpp
OBJ = obj/main.o obj/reader.o obj/compare.o obj/tree.o

sop_trace: $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

obj/%.o: src/%.cpp
	@mkdir -p obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf obj/*.o sop_trace obj