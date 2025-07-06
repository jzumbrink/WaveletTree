GXX = g++ --std=c++20
CC = gcc
EVAL_FLAGS = -O3 -march=native
SPACE_MEASURE_FLAGS = -Iexternal/malloc_count
MALLOC_COUNT_END = -ldl external/malloc_count/malloc_count.o

wave:
	$(CC) -Wall -Wextra -g -c external/malloc_count/malloc_count.c -o external/malloc_count/malloc_count.o
	$(GXX) $(EVAL_FLAGS) -o bin/wave -I ~/include $(SPACE_MEASURE_FLAGS) -Wno-deprecated-declarations -L ~/lib wave.cpp -lsdsl $(MALLOC_COUNT_END)
	./bin/wave