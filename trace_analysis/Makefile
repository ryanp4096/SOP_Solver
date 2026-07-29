
sop_trace: obj/main.o obj/reader.o obj/compare.o
	g++ -pthread -O3 -o sop_trace $^

obj/main.o: src/main.cpp src/main.hpp
	mkdir -p obj
	g++ -O3 -c src/main.cpp -o obj/main.o

obj/reader.o: src/reader.cpp src/reader.hpp
	mkdir -p obj
	g++ -O3 -c src/reader.cpp -o obj/reader.o

obj/compare.o: src/compare.cpp src/compare.hpp
	mkdir -p obj
	g++ -O3 -c src/compare.cpp -o obj/compare.o

clean:
	rm -rf obj/*.o sop_trace obj