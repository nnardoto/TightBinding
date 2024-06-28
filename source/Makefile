	CPPFLAG = -std=c++11
LIBS = -larmadillo
OBJS = BandCalc.o LoadModel.o Parser.o

TightBinding.x: $(OBJS)
	g++ TightBinding.cpp -o TightBinding.x $(OBJS) $(CPPFLAG) $(LIBS)
	mv TightBinding.x ../bin/

BandCalc.o:
	g++ BandCalc.cpp -o BandCalc.o -c $(CPPFLAG) $(LIBS) 

LoadModel.o: Parser.o
	g++ LoadModel.cpp -c -o LoadModel.o $(CPPFLAG) $(LIBS) 

Parser.o:
	g++ Parser.cpp -o Parser.o -c $(CPPFLAG) $(LIBS)

clean:
	rm *.x *.o
