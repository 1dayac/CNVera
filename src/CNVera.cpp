// ASQGParser.cpp : Defines the entry point for the console application.
//

#include "SGUtil.h"
#include "SGVisitors.h"
#include "SGSearch.h"
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <functional>

int VERBOSE = 0;

void logVerbose(std::string message)
{
  if (VERBOSE)
    std::cerr << message << std::endl;
}


typedef std::vector<std::pair<std::string, int>> ScafVector;

class CompareVertexPtr {
public:
  bool operator()(Vertex* t1, Vertex* t2) // Returns true if t1 is earlier than t2
  {
    if (t1->getUntrustyMetric() < t2->getUntrustyMetric()) return true;
    return false;
  }
};


class ComparePQScaffolds {
public:
  bool operator()(Vertex* t1, Vertex* t2) // Returns true if t1 is earlier than t2
  {
    if (t1->getUntrustyMetric() < t2->getUntrustyMetric()) return true;
    return false;
  }
};

typedef std::pair<std::vector<std::string>, int> Path_SGA;

class ComparePathSGA {
public:
  bool operator()(Path_SGA t1, Path_SGA t2) // Returns true if t1 is earlier than t2
  {
    if (t1.second < t2.second) return true;
    return false;
  }
};

struct PositionAndPercentage
{
  int pos;
  double percentage;
  PositionAndPercentage(int _pos, double _percentage)
    : pos(_pos), percentage(_percentage)
  {}
};

std::vector<std::string> &splitnew(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> splitnew(const std::string &s, char delim) {
  std::vector<std::string> elems;
  splitnew(s, delim, elems);
  return elems;
}

//Method to parse Magnolia output and map initial copynumbers to the graph 
void mapCNtoContigs(StringGraph* graph, std::ifstream& inCN)
{
  std::unordered_map<std::string, size_t> CN;
  std::string temp = "";
  while (std::getline(inCN, temp))
  {
    std::string contigName = temp.substr(0, temp.find_first_of("\t"));
    size_t CNEstimation = std::atoi(temp.substr(temp.find_last_of("\t"), temp.length() - temp.find_last_of("\t") + 1).c_str());
    CN[contigName] = CNEstimation;
  }

  auto allVertices = graph->getAllVertices();
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    allVertices[i]->setEstimatedCN(CN[allVertices[i]->getID()]);
    if (allVertices[i]->getSeqLen() > 500 && allVertices[i]->getEstimatedCN() != 0)
    {
      allVertices[i]->setTrustedCN(true);
    }
  }
}

//Storage for distance estimation information
std::unordered_map<std::string, ScafVector> forwardScaf;
std::unordered_map<std::string, ScafVector> reverseScaf;

//Output result
void writeCopyNumbers(StringGraph* graph)
{
  std::ofstream out("new_copynumbers.txt");
  auto allVertices = graph->getAllVertices();
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    out << allVertices[i]->getID() << "\t" << allVertices[i]->getEstimatedCN() << std::endl;
  }
}

//Parse distance estimation file
std::vector<StringVector> parseScaffoldFile(std::ifstream& inScaf)
{
  std::vector<StringVector> scaffolds;
  std::string temp;
  while (std::getline(inScaf, temp))
  {
    StringVector a = split(temp, '\t');
    StringVector res;

    for (int i = 0; i < a.size(); ++i)
    {
      res.push_back(split(a[i], ',')[0]);
    }
    scaffolds.push_back(res);
  }
  return scaffolds;
}


//Parse scaffold file
void parseDEScaffoldFile(std::ifstream& inDEScaf)
{
  std::string temp;
  std::string temp2;
  while (inDEScaf >> temp)
  {
    std::getline(inDEScaf, temp2);

    StringVector a = split(temp2, ';');

    for (int i = 0; i < 2; ++i)
    {
      if (a[i] == "")
        continue;
      StringVector res1 = split(a[i], ' ');
      for (int j = 0; j < res1.size(); ++j)
      {
        if (res1[j] == "")
          continue;
        StringVector res2 = split(res1[j], ',');
        if (i == 0)
        {
          forwardScaf[temp].push_back(std::make_pair(res2[0].substr(0, res2[0].length() - 1), std::atoi(res2[1].c_str())));
        }
        else
        {
          reverseScaf[temp].push_back(std::make_pair(res2[0].substr(0, res2[0].length() - 1), std::atoi(res2[1].c_str())));
        }
      }
    }
  }
}

//Hardcoded value, should work well with any IS < 400
const int insertSize = 400;

std::unordered_map<std::string, VertexPtrVec> reverseAdjList;
std::unordered_map<std::string, VertexPtrVec> forwardAdjList;


struct blastRecord
{
  std::string contigName;
  int start;
  int end;
  bool rc;
  blastRecord(std::string _contigName, int _start, int _end)
    : contigName(_contigName), start(_start), end(_end), rc(false)
  {
    if (start > end)
    {
      rc = true;
      std::swap(start, end);
    }
  }
};

class CompareBlastRec {
public:
  bool operator()(blastRecord b1, blastRecord b2) // Returns true if t1 is earlier than t2
  {
    if (b1.start < b2.start) return true;
    return false;
  }
};

typedef std::vector<blastRecord> blastVector;


void parseBlastFile(blastVector& blastMap, StringGraph* graph)
{
  std::ifstream blastParse("results.txt");
  std::string temp;
  while (getline(blastParse, temp))
  {
    StringVector blastRow = split(temp, '\t');
    if (std::stod(blastRow[2]) > 95.0 && std::stod(blastRow[3].c_str())  > 0.95 * double(graph->getVertex(blastRow[0])->getSeqLen()))
    {
      blastMap.push_back(blastRecord(blastRow[0], std::atoi(blastRow[8].c_str()), std::atoi(blastRow[9].c_str())));
    }
  }
}

//If pathToFind is contained in path, with start in pNode
bool comparepaths(VertexPtrVec pathToFind, pathNode* pNode, bool isForward)
{
  auto currentNode = pNode;
  if ((isForward && pNode->realForward == FORWARD) || (!isForward && pNode->realForward == REVERSE))
  {
    for (int i = 0; i < pathToFind.size(); ++i)
    {
      if (currentNode->current == pathToFind[i])
      {
        if (currentNode->forward != 0)
          currentNode = currentNode->forward;
        else
        {
          if (i != pathToFind.size() - 1)
            break;
        }
      }
      else
      {
        break;
      }
      if (i == pathToFind.size() - 1)
      {
        return true;
      }
    }
  }
  else
  if ((!isForward && pNode->realReverse == REVERSE) || (isForward && pNode->realReverse == FORWARD))
  {
    for (int i = 0; i < pathToFind.size(); ++i)
    {
      if (currentNode->current == pathToFind[i])
      {
        if (currentNode->reverse != 0)
          currentNode = currentNode->reverse;
        else
        {
          if (i != pathToFind.size() - 1)
            break;
        }
      }
      else
      {
        break;
      }
      if (i == pathToFind.size() - 1)
      {
        return true;
      }
    }
  }
  return false;
}


int extendpaths(SGWalkVector& allpaths, bool isForward)
{
  int totalPathSize = 0;
  for (int i = 0; i < allpaths.size(); ++i)
  {

    auto path = allpaths[i];
    auto vertexSeq = path.getVertices();
    auto firstVertex = vertexSeq[0];
    int currpathsSize = 0;
    for (auto d : firstVertex->getPaths())
    {
      if (comparepaths(vertexSeq, d, isForward))
      {
        currpathsSize++;
      }
    }
    totalPathSize += std::max(1, currpathsSize);
  }
  return totalPathSize;
}

enum Contained
{
  FIRSTINSECOND, SECONDINFIRST, NONE
};

//Which of path is contained in another?
Contained isContained(SGWalk walk1, SGWalk walk2)
{
  Contained answer = SECONDINFIRST;
  if (walk1.getEndToStartDistance() < walk2.getEndToStartDistance())
  {
    answer = FIRSTINSECOND;
    std::swap(walk1, walk2);
  }

  for (int i = 0; i < std::min(walk1.getEdgeIndex().size(), walk2.getEdgeIndex().size()); ++i)
  {
    if (walk1.getEdgeIndex()[i] != walk2.getEdgeIndex()[i])
    {
      return NONE;
    }
  }
  return answer;
}

std::unordered_map<std::string, int> CNbyNeighbours;


//Remove uncontained paths and then call split extend algorithm
int uncontainedpaths_SGA(Vertex* startID, ScafVector& scaffolds, bool isForward, StringGraph* graph)
{
  int answer = 0;
  SGWalkVector superWalks;
  for (int i = 0; i < scaffolds.size(); ++i)
  {
    SGWalkVector outWalks;
    SGSearch::findWalks(startID, graph->getVertex(scaffolds[i].first), isForward ? ED_SENSE : ED_ANTISENSE, insertSize, 3000, false, outWalks);
    logVerbose("Number of paths from " + startID->getID() + " to " + scaffolds[i].first + " is " + std::to_string(outWalks.size()));
    int difference = 100000;
    int bestPathIndex = 0;
    if (outWalks.size() != 0)
    {
      SGWalkVector setOfContenders;

      for (int j = 0; j < outWalks.size(); ++j)
      {
        if (abs(outWalks[j].getEndToStartDistance() - scaffolds[i].second - 12) < difference)
        {
          difference = abs(outWalks[j].getEndToStartDistance() - scaffolds[i].second - 12);
          bestPathIndex = j;
        }
        if (abs(outWalks[j].getEndToStartDistance() - scaffolds[i].second - 12) < 100)
        {
          for (auto d : startID->getPaths())
          {
            if (comparepaths(outWalks[j].getVertices(), d, isForward))
            {
              setOfContenders.push_back(outWalks[j]);
              break;
            }
          }
        }
      }

      if (setOfContenders.size() == 0 && difference < 100)
      {
        logVerbose("Don't find any contender. Index of best path - " + std::to_string(bestPathIndex));
        superWalks.push_back(outWalks[bestPathIndex]);
      }
      else
      {
        logVerbose("Number of contenders for this pair - " + std::to_string(setOfContenders.size()));
        for (auto d : setOfContenders)
        {
          superWalks.push_back(d);
        }
      }
    }
  }

  std::vector<bool> checked(superWalks.size());
  for (int i = 0; i < superWalks.size(); ++i)
  {

    for (int j = i + 1; j < superWalks.size(); ++j)
    {
      if (checked[i] || checked[j])
        continue;
      switch (isContained(superWalks[i], superWalks[j]))
      {
        case NONE:
          break;
        case FIRSTINSECOND:
          checked[i] = true;
          break;
        case SECONDINFIRST:
          checked[j] = true;
          break;
      }
    }
  }

  std::unordered_map<std::string, int> CNByNeighboursTemp;
  for (size_t i = 0; i < superWalks.size(); ++i)
  {
    if (!checked[i])
    {
      for (size_t j = 0; j < superWalks[i].getVertices().size(); ++j)
      {
        CNByNeighboursTemp[superWalks[i].getVertices()[j]->getID()]++;
      }
    }
  }


  SGWalkVector allpaths;

  for (size_t i = 0; i < superWalks.size(); ++i)
  {
    if (!checked[i])
    {
      answer++;
      allpaths.push_back(superWalks[i]);
    }
  }

  logVerbose("Total number of uncontained paths - " + std::to_string(answer));

  for (auto d : CNByNeighboursTemp)
  {
    if (d.second > CNbyNeighbours[d.first])
      CNbyNeighbours[d.first] = d.second;
  }

  int res = extendpaths(allpaths, isForward);
  logVerbose("Number of paths after extension step - " + std::to_string(res));

  return res;
}


//Map reference to graph
void mappathsToGraph(StringGraph * graph, blastVector& blastVector)
{
  std::sort(blastVector.begin(), blastVector.end(), CompareBlastRec());
  pathNode* start = new pathNode(graph->getVertex(blastVector[0].contigName));
  graph->getVertex(blastVector[0].contigName)->addPath(start);

  for (int i = 0; i < blastVector.size() - 1; ++i)
  {
    if (blastVector[i].end - blastVector[i + 1].start >= 55)
    {
      pathNode * curr = new pathNode(graph->getVertex(blastVector[i + 1].contigName), graph->getVertex(blastVector[i].contigName)->getPaths()[graph->getVertex(blastVector[i].contigName)->getPaths().size() - 1]);
      if (blastVector[i + 1].rc == true)
      {
        curr->realReverse = FORWARD;
        curr->realForward = REVERSE;
      }
      else
      {
        curr->realReverse = REVERSE;
        curr->realForward = FORWARD;
      }
      graph->getVertex(blastVector[i + 1].contigName)->addPath(curr);
      graph->getVertex(blastVector[i].contigName)->setForwardToPath(graph->getVertex(blastVector[i + 1].contigName)->getPaths()[graph->getVertex(blastVector[i + 1].contigName)->getPaths().size() - 1]);
    }
    else
    {
      pathNode * curr = new pathNode(graph->getVertex(blastVector[i + 1].contigName));
      graph->getVertex(blastVector[i + 1].contigName)->addPath(curr);
    }
  }
}



int main_work(int argc, char* argv[])
{
  // Arguments
  // primary-graph.asqg -- SGA graph structure
  // copynumbers.txt -- Magnolia solution
  // libPE.de -- contig supported by paired reads information
  // scaffold.scaf -- scaffolding information
  // reference.fa -- reference
  // assemblem100-contigs.fa Contigs file

  std::string filename = argv[1];
  std::string referenceFilename = argv[5];
  std::string contigsFilename = argv[6];
  std::string makeblastdbCommand = "$SHELL -c 'makeblastdb -dbtype nucl -in " + referenceFilename + " -out tempDatabase'";
  //make blast database for reference
  int excode = system(makeblastdbCommand.c_str());
  std::string blastContigsCommand = "$SHELL -c 'blastn -db tempDatabase -query " + contigsFilename + " -outfmt 6 -dust yes -word_size 80 -evalue 10 -perc_identity 95 -out results.txt'";
  system(blastContigsCommand.c_str());

  if (excode)
  {
    std::cerr << "Error making blast db" << std::endl;
    return 1;
  }

  //contruct adjacency lists
  StringGraph* graph = SGUtil::loadASQG(filename, 0);
  blastVector blastOfAllContigs;
  parseBlastFile(blastOfAllContigs, graph);
  mappathsToGraph(graph, blastOfAllContigs);



  std::ifstream CNEstimation(argv[2]);

  auto allVertices = graph->getAllVertices();

  std::ifstream scafParser(argv[3]);

  parseDEScaffoldFile(scafParser);

  std::ifstream scaffoldFile(argv[4]);

  mapCNtoContigs(graph, CNEstimation);

  for (auto d : graph->getAllVertices())
  {
    std::cout << d->getID() << " is processed" << std::endl;
    std::string idname = d->getID();
    if (d->getSeqLen() > 10000)
    {
      continue;
    }

    int forwardCount = 1;
    int reverseCount = 1;
    bool isTrulyProcessed = false;
    if (forwardScaf.find(d->getID()) != forwardScaf.end())
    {

      forwardCount = uncontainedpaths_SGA(graph->getVertex(d->getID()), forwardScaf[d->getID()], true, graph);
      logVerbose(d->getID() + " forward count is " + std::to_string(forwardCount));

      if (forwardCount)
        isTrulyProcessed = true;
    }
    if (reverseScaf.find(d->getID()) != reverseScaf.end())
    {
      reverseCount = uncontainedpaths_SGA(graph->getVertex(d->getID()), reverseScaf[d->getID()], false, graph);
      logVerbose(d->getID() + " reverse count is " + std::to_string(reverseCount));
      if (reverseCount)
      {
        isTrulyProcessed = true;
      }
    }

    if (isTrulyProcessed)
    {
      d->setEstimatedWithPairedReads(true);
    }

    if (std::max(forwardCount, reverseCount) > d->getEstimatedCN())
    {
      d->setEstimatedCN(std::max(std::max(forwardCount, reverseCount), 1));
      logVerbose(d->getID() + " is changed on " + std::to_string(std::max(std::max(forwardCount, reverseCount), 1)));
       d->setTrustedCN(false);
    }
    else
    {
      logVerbose(d->getID() + " is not changed");
    }
  }

  for (auto d : graph->getAllVertices())
  {
    if (d->getSeqLen() > 800)
    {
      continue;
    }

    if (!d->getEstimatedWithPairedReads() && d->getEstimatedCN() < CNbyNeighbours[d->getID()])
    {
      logVerbose(d->getID() + " is changed on " + std::to_string(CNbyNeighbours[d->getID()]) + " by neighbours");
      d->setEstimatedCN(CNbyNeighbours[d->getID()]);
    }
  }

  // Write the results
  writeCopyNumbers(graph);
  SGFastaVisitor av("new_contigs.fasta");
  graph->visit(av);
  return 0;
}

int main(int argc, char* argv[])
{
  main_work(argc, argv);
}

