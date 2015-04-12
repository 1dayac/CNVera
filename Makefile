CXX=g++
CXXFLAGS= -std=c++11 -O3 -Wl,--stack=256000000

all:
	$(CXX) $(CXXFLAGS) Alphabet.cpp ASQG.cpp ASQG_main.cpp Bigraph.cpp BitChar.cpp DNAString.cpp doppelGraph.cpp Edge.cpp EdgeDesc.cpp  Interval.cpp Match.cpp MultiOverlap.cpp Quality.cpp QualityVector.cpp SeqCoord.cpp SGAlgorithms.cpp SGSearch.cpp SGUtil.cpp SGWalk.cpp SGVisitors.cpp SQG.cpp Util.cpp Vertex.cpp -o ASQGParser

Vertex.o: Vertex.cpp
	$(CXX) $(CXXFLAGS) -c Vertex.cpp

Utils.o: Utils.cpp
	$(CXX) $(CXXFLAGS) -c Utils.cpp

.PHONY: clean
clean:
	rm -f ASQGParser
