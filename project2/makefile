CPPflags = g++ -c -std=c++11
LIB =  -larmadillo
PROG = main

all: compile link clean

compile:
	${CPPflags} -c ${PROG}.cpp

link:
	g++ ${PROG}.o -o ${PROG}.exe ${LIB}

run:
	./${PROG}.exe

plot:
	python3 ./plotfiles/plot6.py
	python3 ./plotfiles/plot7.py

clean:
	rm -f *.o *~
