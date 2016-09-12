binPath=./bin
cpp=./src/ppso.cpp
all:compile run

compile:$(cpp)
	mpic++.mpich -O2 -o $(binPath)/LPSO -DALGORITHM=0 $(cpp)
	mpic++.mpich -O2 -o $(binPath)/GPSO -DALGORITHM=1 $(cpp)
	mpic++.mpich -O2 -o $(binPath)/BPSO -DALGORITHM=2 $(cpp)
	mpic++.mpich -O2 -o $(binPath)/CLPSO -DALGORITHM=3 $(cpp)
	mpic++.mpich -O2 -o $(binPath)/PPSO -DALGORITHM=4 $(cpp)

run:$(binPath)/*
	cd $(binPath) && make
	@echo ****Result is saved in ./output/all.txt******



