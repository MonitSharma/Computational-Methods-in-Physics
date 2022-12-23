# use this by doing this in terminal:
# make
# make run (or make compile, make clean etc.)
CPPflags = g++ -std=c++11 -O2
LIB =  -larmadillo
all: compile link clean

compile:
	${CPPflags} -c main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/#${LIB}

link:
	${CPPflags} -o main.exe main.o PenningTrap.o Particle.o ${LIB}

run:
	./main.exe

clean:
	rm -f *.o *~

plot:
	python3 ./plotfiles/9_1.py
	python3 ./plotfiles/9_2.py 100
	python3 ./plotfiles/9_2.py 1000
	python3 ./plotfiles/9_3.py 100
	python3 ./plotfiles/9_3.py 1000
	python3 ./plotfiles/9_4.py 100
	python3 ./plotfiles/9_4.py 1000
	python3 ./plotfiles/9_5.py
	python3 ./plotfiles/10_6.py
	#python3 ./plotfiles/10_7.py
	python3 ./plotfiles/10_7.py average
