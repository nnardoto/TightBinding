unamed.x:
	g++ unamed.cpp -o unamed.x BandCalc.o -std=c++11 -larmadillo
BandCalc.o:
	g++ BandCalc.cpp -o BandCalc.o -c -std=c++11 -larmadillo
	 
