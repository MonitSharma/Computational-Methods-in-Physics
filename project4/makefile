# use this by doing this in terminal:
# make
# make run (or make compile, make clean etc.)
CPPflags = g++ -std=c++11 -O2  # Mac users add -Xpreprocessor before -fopenmp
#LIB =  -larmadillo
all: compile link clean

compile:
	${CPPflags} -c main.cpp src/Ising.cpp -Xpreprocessor -fopenmp -I include/ # ${LIB}

link:
	${CPPflags} -o main.exe main.o Ising.o -larmadillo -lomp # ${LIB}

run:
	./main.exe

clean:
	rm -f *.o *~

plot:
	python3 ./plotfiles/problem4.py
	python3 ./plotfiles/problem5.py
	#python3 ./plotfiles/problem6.py
	python3 ./plotfiles/problem8.py
	python3 ./plotfiles/problem9.py



# CPPflags = c++ #-O3 -std=c++11
# LIB =  #-larmadillo # -llapack -lblas
# PROG = main
#
# all: compile link clean #run
#
# compile:
# 	${CPPflags} -c ${PROG}.cpp functions.cpp ${LIB}
#
# link:
# 	${CPPflags} -o ${PROG}.out ${PROG}.o functions.o ${LIB}
#
# run:
# 	./${PROG}.out
#
# clean:
# 	rm -f *.o *~





#
#
#
# # --- Alternative: wildcard including all ---#
# all_alternative: compile_alternative link_alternative clean
#
# compile_alternative:
# 	${CPPflags} -c ${PROG}.cpp $(wildcard *.cpp) ${LIB}
#
# link_alternative:
# 	${CPPflags} -o ${PROG}.out $(wildcard *.o) ${LIB}
#
#
# # --- Compact: more simple ---#
# all_compact: compile_compact
#
# compile_compact:
# 	c++ -o main.out main.cpp functions.cpp
#


#--- Old method ---#
# CPPflags = c++ #-O3 -std=c++11
# LIB =  #-larmadillo -llapack -lblas
# PROG = main
#
#
# ${PROG}: 						${PROG}.o functions.o
# 										${CPPflags} ${PROG}.o functions.o ${LIB} -o ${PROG}.exe
#
#
# ${PROG}.o: 					${PROG}.cpp
# 										${CPPflags} -c ${PROG}.cpp
#
# functions.o: 				functions.cpp functions.hpp
# 										${CPPflags} -c functions.cpp
#
#
# clean:
# 	rm -f *.o *~
