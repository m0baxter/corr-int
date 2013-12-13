
COMP = g++
FLAGS = -std=c++11 -Wall -pedantic -flto -Ofast -funroll-loops -march=native
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
	rm  *.d *.o
	
cleanall:
	rm  *.d *.o main
	
-include $(DPND)
