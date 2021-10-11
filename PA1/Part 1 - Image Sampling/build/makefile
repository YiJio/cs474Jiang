#Compiler Flags
CC = g++
CXXFLAGS = -pedantic -std=c++14 -Wall -g

# .o Compilation
O_FILES=main.o image.o ReadImage.o ReadImageHeader.o WriteImage.o ImageSampling.o

# Point to includes of local directories
INDLUDES=-I../include


all: $(O_FILES)
	$(CC) $(CXXFLAGS) -o main $(O_FILES) $(LIBS) #$(PATHI) $(PATHL)

main.o: ../src/main.cpp
	$(CC) $(CXXFLAGS) -c ../src/main.cpp -o main.o $(INDLUDES)

image.o: ../src/image.cpp
	$(CC) $(CXXFLAGS) -c ../src/image.cpp -o image.o $(INDLUDES)

ReadImage.o: ../src/ReadImage.cpp
	$(CC) $(CXXFLAGS) -c ../src/ReadImage.cpp -o ReadImage.o $(INDLUDES)

ReadImageHeader.o: ../src/ReadImageHeader.cpp
	$(CC) $(CXXFLAGS) -c ../src/ReadImageHeader.cpp -o ReadImageHeader.o $(INDLUDES)

WriteImage.o: ../src/WriteImage.cpp
	$(CC) $(CXXFLAGS) -c ../src/WriteImage.cpp -o WriteImage.o $(INDLUDES)

ImageSampling.o: ../src/ImageSampling.cpp
	$(CC) $(CXXFLAGS) -c ../src/ImageSampling.cpp -o ImageSampling.o $(INDLUDES)

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
