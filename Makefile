CC = clang++
CSTD = -std=c++14

debug:
	$(CC) -g -o main src/main.cpp

all: main

main: src/main.cpp include/*.h
	$(CC) $(CSTD) src/main.cpp -o main

clean:
	rm -f main