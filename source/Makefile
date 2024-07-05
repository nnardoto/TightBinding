CPPFLAG = -std=c++23
LIBS = -larmadillo
OBJS = BandCalc.o LoadModel.o Parser.o Export2Xatu.o KPath.o

TightBinding.x: $(OBJS)
	g++ TightBinding.cpp -o TightBinding.x $(OBJS) $(CPPFLAG) $(LIBS)

BandCalc.o:
	g++ BandCalc.cpp -o BandCalc.o -c $(CPPFLAG) $(LIBS) 

LoadModel.o: Parser.o
	g++ LoadModel.cpp -c -o LoadModel.o $(CPPFLAG) $(LIBS) 

Parser.o:
	g++ Parser.cpp -o Parser.o -c $(CPPFLAG) $(LIBS)

Export2Xatu.o:
	g++ Export2Xatu.cpp -o Export2Xatu.o -c $(CPPFLAG) $(LIBS)

KPath.o: BandCalc.o
	g++ KPath.cpp -o KPath.o -c $(CPPFLAG) $(LIBS)

clean:
	rm *.x *.o
