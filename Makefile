GXX = g++ --std=c++20
EVAL_FLAGS = -O3 -march=native

wave:
	$(GXX) $(EVAL_FLAGS) -o bin/wave -I ~/include -Wno-deprecated-declarations -L ~/lib wave.cpp -lsdsl
	./bin/wave

sdsl-wave:
	$(GXX) $(EVAL_FLAGS) -o bin/sdsl-wave -I ~/include -Wno-deprecated-declarations -L ~/lib sdsl-wave.cpp -lsdsl