
COMP = g++
FLAGS = -Wall -pedantic

main : main.o wavefunction.o AngularMomentum.o newtoncotes.o makefile
	      $(COMP) -o main main.o wavefunction.o AngularMomentum.o newtoncotes.o
	      
newtoncotes.o : newtoncotes.cpp newtoncotes.hpp makefile
		          $(COMP) -c $(FLAGS) newtoncotes.cpp

AngularMomentum.o : AngularMomentum.cpp AngularMomentum.hpp makefile
				    $(COMP) -c $(FLAGS) AngularMomentum.cpp

wavefunction.o : wavefunction.cpp wavefunction.hpp AngularMomentum.hpp newtoncotes.hpp makefile
	                $(COMP) -c $(FLAGS) wavefunction.cpp

main.o : main.cpp wavefunction.hpp makefile
	        $(COMP) -c $(FLAGS) main.cpp
clean:
		rm  main.o wavefunction.o AngularMomentum.o newtoncotes.o
