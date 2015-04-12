CXX=g++
CXXFLAGS= -std=c++11 -O3
PATHTOFILES=./src
PATHTOFILESSIM=./tools/SVSim

all:
	$(CXX) $(CXXFLAGS) $(PATHTOFILES)/Alphabet.cpp $(PATHTOFILES)/ASQG.cpp $(PATHTOFILES)/CNVera.cpp $(PATHTOFILES)/Bigraph.cpp $(PATHTOFILES)/BitChar.cpp $(PATHTOFILES)/DNAString.cpp $(PATHTOFILES)/Edge.cpp $(PATHTOFILES)/EdgeDesc.cpp  $(PATHTOFILES)/Interval.cpp $(PATHTOFILES)/Match.cpp $(PATHTOFILES)/MultiOverlap.cpp $(PATHTOFILES)/Quality.cpp $(PATHTOFILES)/QualityVector.cpp $(PATHTOFILES)/SeqCoord.cpp $(PATHTOFILES)/SGAlgorithms.cpp $(PATHTOFILES)/SGSearch.cpp $(PATHTOFILES)/SGUtil.cpp $(PATHTOFILES)/SGWalk.cpp $(PATHTOFILES)/SGVisitors.cpp $(PATHTOFILES)/SQG.cpp $(PATHTOFILES)/Util.cpp $(PATHTOFILES)/Vertex.cpp -o CNVera

  $(CXX) $(CXXFLAGS) $(PATHTOFILESSIM)/fasta.cpp $(PATHTOFILESSIM)/SVSimulator.cpp -o $(PATHTOFILESSIM)/SVSim
Vertex.o: Vertex.cpp
	$(CXX) $(CXXFLAGS) -c Vertex.cpp

Utils.o: Utils.cpp
	$(CXX) $(CXXFLAGS) -c Utils.cpp

.PHONY: clean
clean:
	rm -f CNVera
	rm -f CNVera.exe
