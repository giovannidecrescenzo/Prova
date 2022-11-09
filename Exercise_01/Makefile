CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o functions.o analysis.o
	$(CC) random.o main.o functions.o analysis.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : ../Random/random.cpp ../Random/random.h
	$(CC) -c ../Random/random.cpp -o random.o $(CFLAGS)
functions.o : ../Functions/functions.cpp ../Functions/functions.h
	$(CC) -c ../Functions/functions.cpp -o functions.o $(CFLAGS)
analysis.o : ../Functions/analysis.cpp ../Functions/analysis.h
	$(CC) -c ../Functions/analysis.cpp -o analysis.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.txt
exec :
	./main.exe
