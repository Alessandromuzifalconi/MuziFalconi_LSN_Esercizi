CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistics.o individual.o 
	$(CC) statistics.o random.o individual.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
individual.o : individual.cpp individual.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
