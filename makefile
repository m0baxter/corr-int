
COMP = g++
FLAGS = -Wall -pedantic -O3 -march=native
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
DPND = $(SRC:.cpp=.d)

main : $(OBJ)
	$(COMP) $(FLAGS) $^ -o $@
	
%.o: %.cpp
	$(COMP) $(FLAGS) -MMD -MP -c $< -o $@

main.o : main.cpp wavefunction.hpp makefile
	        $(COMP) -c $(FLAGS) main.cpp
clean:
	rm  *.o main
	
-include $(DPND)
