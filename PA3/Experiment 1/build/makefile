#Compiler Flags
CC = g++
CXXFLAGS = -pedantic -std=c++14 -Wall -g

# .o Compilation
O_FILES=main.o fft.o

# Point to includes of local directories
INDLUDES=-I../include


all: $(O_FILES)
	$(CC) $(CXXFLAGS) -o main $(O_FILES) $(LIBS) #$(PATHI) $(PATHL)

main.o: ../src/main.cpp
	$(CC) $(CXXFLAGS) -c ../src/main.cpp -o main.o $(INDLUDES)

fft.o: ../src/fft.cpp
	$(CC) $(CXXFLAGS) -c ../src/fft.cpp -o fft.o $(INDLUDES)

clean:
	rm -f main
	rm -f *.o
	reset
	@echo
	@echo
	@echo
	@echo Make clean, removed object and executable files ...
	@echo Contents in this Directory ...
	ls -lt
	@echo
	@echo
	@echo
