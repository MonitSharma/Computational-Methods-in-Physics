CPPflags = g++ -c -std=c++11
LIB =  -larmadillo
PROG = main

all: compile link clean #run

compile:
	${CPPflags} -c ${PROG}.cpp

link:
	g++ ${PROG}.o -o ${PROG}.exe ${LIB}

run:
	# main.exe {N} where N is number of steps
	./${PROG}.exe 10
	./${PROG}.exe 100
	./${PROG}.exe 1000
	./${PROG}.exe 10000
	./${PROG}.exe 100000
	./${PROG}.exe 1000000

plot:
	python3 ./plotfiles/plot_exact.py
	python3 ./plotfiles/general_vs_exact.py
	python3 ./plotfiles/errorplot.py

clean:
	rm -f *.o *~
