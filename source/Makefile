CPPFLAG = -std=c++23
LIBS = -larmadillo
OBJS = TBModel.o BandCalc.o LoadModel.o Parser.o Export2Xatu.o KPath.o Symmetrize.o

TightBinding.x: TightBinding.cpp $(OBJS)
	g++ TightBinding.cpp -o TightBinding.x $(OBJS) $(CPPFLAG) $(LIBS)

TBModel.o: TBModel.h TBModel.cpp
	g++ TBModel.cpp -o TBModel.o -c $(CPPFLAG) $(LIBS)  

BandCalc.o: BandCalc.cpp TBModel.o
	g++ BandCalc.cpp -o BandCalc.o -c $(CPPFLAG) $(LIBS) 

LoadModel.o: LoadModel.cpp Parser.o TBModel.o
	g++ LoadModel.cpp -c -o LoadModel.o $(CPPFLAG) $(LIBS) 

Parser.o: Parser.cpp TBModel.o
	g++ Parser.cpp -o Parser.o -c $(CPPFLAG) $(LIBS)

Export2Xatu.o: Export2Xatu.cpp TBModel.o
	g++ Export2Xatu.cpp -o Export2Xatu.o -c $(CPPFLAG) $(LIBS)

KPath.o: KPath.cpp BandCalc.o TBModel.o
	g++ KPath.cpp -o KPath.o -c $(CPPFLAG) $(LIBS)

Symmetrize.o: Symmetrize.cpp TBModel.o
	g++ Symmetrize.cpp -o Symmetrize.o -c $(CPPFLAG) $(LIBS)

clean:
	rm *.x *.o
