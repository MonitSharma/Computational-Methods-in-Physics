# use this by doing this in terminal:
# make
# make run (or make compile, make clean etc.)

CPPflags = g++ -std=c++14 -O2
LIB =  -larmadillo
all: compile link clean

compile:
	${CPPflags} -c main.cpp src/Crank.cpp -Xpreprocessor -fopenmp -I include/ # ${LIB}

link:
	${CPPflags} -o main.exe main.o Crank.o -larmadillo -lomp #${LIB}

run:
	./main.exe "./inputs/input.txt"
extra:
	./main.exe "./inputs/extra.txt"
clean:
	rm -f *.o *~

plot:
	#Animation format is: animation.py 'filename' 'heatmap label' 'number of slits in overlay' 'number of snapshots' 'snapshot times'
	python3 plotfiles/prob_vs_t.py 'Problem_7.10_outputCube_slits_0.dat' 0
	python3 plotfiles/prob_vs_t.py 'Problem_7.20_outputCube_slits_2.dat' 2
	python3 plotfiles/animation.py 'Problem_7.10_outputCube_slits_0.dat' 'p(x,y,t)' 0 0
	python3 plotfiles/animation.py 'Problem_7.20_outputCube_slits_2.dat' 'p(x,y,t)' 2 0
	python3 plotfiles/animation.py 'Problem_8.00_outputCube_slits_2.dat' 'p(x,y,t)' 2 3 0.0 0.001 0.002
	python3 plotfiles/animation.py 'Problem_8.00_ReIm_outputCube_slits_2.dat' 'Real(u)' 2 3 0.0 0.001 0.002
	python3 plotfiles/animation.py 'Problem_8.00_ReIm_outputCube_slits_2.dat' 'Imag(u)' 2 3 0.0 0.001 0.002
	python3 plotfiles/detection_prob.py 'Problem_9.10_outputMat_slits_2.dat' 2
	python3 plotfiles/detection_prob.py 'Problem_9.21_outputMat_slits_1.dat' 1
	python3 plotfiles/detection_prob.py 'Problem_9.22_outputMat_slits_3.dat' 3



	# python3 plotfiles/animation.py 'Problem_8.00_outputCube_slits_2.dat' 0.00 0.5 1.00
