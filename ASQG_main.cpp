// ASQGParser.cpp : Defines the entry point for the console application.
//

#include "SGUtil.h"
#include "doppelGraph.h"
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

int z = 0;
int tryToDelete = 0;

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

void deleteEdges(std::unordered_map<std::string, std::vector<PositionAndPercentage>>& blastRes, Vertex* v, std::unordered_map<std::string, VertexPtrVec>& reverseAdjList, std::unordered_map<std::string, VertexPtrVec>& forwardAdjList)
{
  if (blastRes[v->getID()].size() == 0)
    return;
  int positionOfUnique = blastRes[v->getID()][0].pos;
  auto edges = forwardAdjList[v->getID()];
  auto revEdges = reverseAdjList[v->getID()];
  std::unordered_set<std::string> goodFwdID; 
  bool isFwdGood = false;
  double bestPercentage = 0.0;
  int bestContigLength = 0;

  if (edges.size() >= 2)
  {
    for (size_t i = 0; i < edges.size(); ++i)
    {
      for (size_t j = 0; j < blastRes[edges[i]->getID()].size(); ++j)
      {
        if (std::abs((int)blastRes[edges[i]->getID()][j].pos - (int)positionOfUnique) < 200 || std::abs((int)blastRes[edges[i]->getID()][j].pos - (int)positionOfUnique + (int)v->getSeqLen()) < 200)
        {
          isFwdGood = true;

          //some bad things
          if (bestPercentage < blastRes[edges[i]->getID()][j].percentage && bestContigLength < edges[i]->getSeqLen() * 5)
          {
            bestContigLength = edges[i]->getSeqLen();
            bestPercentage = blastRes[edges[i]->getID()][j].percentage;
            goodFwdID.clear();
            goodFwdID.insert(edges[i]->getID());
          }
        }
      }
    }

  }
  else
  {
    if (edges.size() == 1)
      goodFwdID.insert(edges[0]->getID());
  }

  bool isRevGood = false;
  std::unordered_set<std::string> goodRevID;
  bestPercentage = 0.0;
  bestContigLength = 0;

  if (revEdges.size() >= 2)
  {

    for (int i = 0; i < revEdges.size(); ++i)
    {
      for (int j = 0; j < blastRes[revEdges[i]->getID()].size(); ++j)
      {
        if (std::abs((int)blastRes[revEdges[i]->getID()][j].pos - (int)positionOfUnique + (int)v->getSeqLen()) < 200 || std::abs((int)blastRes[revEdges[i]->getID()][j].pos - (int)positionOfUnique) < 200)
        {
          isRevGood = true;
          
          if (bestPercentage < blastRes[revEdges[i]->getID()][j].percentage && bestContigLength < revEdges[i]->getSeqLen() * 5)
          {
            bestContigLength = revEdges[i]->getSeqLen();
            bestPercentage = blastRes[revEdges[i]->getID()][j].percentage;
            goodRevID.clear();
            goodRevID.insert(revEdges[i]->getID());
          }
        }
      }
    }
  }
  else
  {
    if (revEdges.size() == 1)
      goodRevID.insert(revEdges[0]->getID());
  }

  if (isFwdGood || isRevGood)
  {
    for (int i = 0; i < v->getEdges().size(); ++i)
    {
      if (goodFwdID.find(v->getEdges()[i]->getEndID()) == goodFwdID.end() && goodRevID.find(v->getEdges()[i]->getEndID()) == goodRevID.end())
      {
        z++;
        if (v->getEdges()[i]->getTwin() != NULL)
           v->getEdges()[i]->getEnd()->deleteEdge(v->getEdges()[i]->getTwin());
        v->deleteEdge(v->getEdges()[i]);
      }
    }
  }
}



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

  std::ofstream outTrusted("trusted.txt");
  std::ofstream outUntrusted("untrusted.txt");

  auto allVertices = graph->getAllVertices();
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    allVertices[i]->setEstimatedCN(CN[allVertices[i]->getID()]);
    if (allVertices[i]->getSeqLen() > 500 && allVertices[i]->getEstimatedCN() != 0)
    {
      allVertices[i]->setTrustedCN(true);
      outTrusted << allVertices[i]->getID() << std::endl;
    }
    else
    {
      outUntrusted << allVertices[i]->getID() << std::endl;
    }   
  }
}

const int overlapAssembly = 110;

std::unordered_map<std::string, ScafVector> forwardScaf;
std::unordered_map<std::string, ScafVector> reverseScaf;

void solveSimpleImplementation()
{

}

void solveSimpleCases(StringGraph* graph)
{





  std::unordered_map<std::string, VertexPtrVec> reverseAdjList;
  std::unordered_map<std::string, VertexPtrVec> forwardAdjList;

  //contruct adjacency lists
  auto allVertices = graph->getAllVertices();
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    SGWalkVector forwardW;
    SGWalkVector reverseW;

    SGSearch::findAllWalks(allVertices[i], ED_SENSE, 1000, 1, forwardW);
    SGSearch::findAllWalks(allVertices[i], ED_ANTISENSE, 1000, 1, reverseW);

    for (size_t j = 0; j < forwardW.size(); ++j)
    {
      forwardAdjList[allVertices[i]->getID()].push_back(forwardW[j].getLastVertex());
    }
    for (size_t j = 0; j < reverseW.size(); ++j)
    {
      reverseAdjList[allVertices[i]->getID()].push_back(reverseW[j].getLastVertex());
    }
  }


  //if subgraph has directed star representation - we could use some very basic methods
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    auto forward = forwardAdjList[allVertices[i]->getID()];
    auto reverse = reverseAdjList[allVertices[i]->getID()];

    if (forward.size() == 0 || reverse.size() == 0)
      continue;

    SGWalkVector forward1;
    SGWalkVector forward2;

    SGWalkVector reverse1;
    SGWalkVector reverse2;

    bool isStar = true;
    for (int j = 0; j < forward.size(); ++j)
    {
      SGSearch::findAllWalks(forward[j], ED_SENSE, 1000, 1, forward1);
      SGSearch::findAllWalks(forward[j], ED_ANTISENSE, 1000, 1, forward2);
      if (!((forward1.size() == 1 && forward1[0].getLastVertex() == allVertices[i]) || (forward2.size() == 1 && forward2[0].getLastVertex() == allVertices[i])))
      {
        isStar = false;
        break;
      }
      forward1.clear();
      forward2.clear();
    }
    if (!isStar)
      continue;

    for (size_t j = 0; j < reverse.size(); ++j)
    {
      SGSearch::findAllWalks(reverse[j], ED_SENSE, 1000, 1, reverse1);
      SGSearch::findAllWalks(reverse[j], ED_ANTISENSE, 1000, 1, reverse2);
      if (!((reverse1.size() == 1 && reverse1[0].getLastVertex() == allVertices[i]) || (reverse2.size() == 1 && reverse2[0].getLastVertex() == allVertices[i])))
      {
        isStar = false;
        break;
      }
      reverse1.clear();
      reverse2.clear();

    }
    if (!isStar)
      continue;
   
    bool amITrusted = allVertices[i]->getTrustedCN();
    
    int numberOfUntrustedForward = 0;
    int numberOfUntrustedReverse = 0;


    int forwardCoverage = 0;
    int backwardCoverage = 0;

    for (auto d : forward)
    {
      forwardCoverage += d->getEstimatedCN();
      if (!d->getTrustedCN())
      {
        numberOfUntrustedForward++;
      }
    }

    for (auto d : reverse)
    {
      backwardCoverage += d->getEstimatedCN();
      if (!d->getTrustedCN())
      {
        numberOfUntrustedReverse++;
      }
    }



    if (forwardCoverage == backwardCoverage && forwardCoverage != allVertices[i]->getEstimatedCN())
    {
      allVertices[i]->setEstimatedCN(forwardCoverage);
      solveSimpleImplementation();
    }
    //else
    //if (forwardCoverage != backwardCoverage && allVertices[i]->getEstimatedCN() < forwardCoverage && allVertices[i]->getEstimatedCN() < backwardCoverage && allVertices[i]->getSeqLen() < 800)
    //{
    //  allVertices[i]->setEstimatedCN(std::min(forwardCoverage, backwardCoverage));
    //}
  }
}

void writeCopyNumbers(StringGraph* graph)
{
  std::ofstream out("new_copynumbers.txt");
  auto allVertices = graph->getAllVertices();
  for (size_t i = 0; i < allVertices.size(); ++i)
  {
    out << allVertices[i]->getID() << "\t" << allVertices[i]->getEstimatedCN() << std::endl;
  }
}


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

const int insertSize = 400;
std::unordered_map<std::string, VertexPtrVec> reverseAdjList;
std::unordered_map<std::string, VertexPtrVec> forwardAdjList;


typedef std::pair<std::vector<std::string>, int> Path_SGA;

class ComparePathSGA {
public:
  bool operator()(Path_SGA t1, Path_SGA t2) // Returns true if t1 is earlier than t2
  {
    if (t1.second < t2.second) return true;
    return false;
  }
};

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


void parseBlastFile(std::unordered_map<int, blastVector>& blastMap)
{
  std::ifstream blastParse("tempResults.txt");
  std::string temp;
  while (getline(blastParse, temp))
  {
    StringVector blastRow = split(temp, '\t');
    StringVector contigInfo = split(blastRow[0], '|');
    if ((blastRow[2] == "100.00" && blastRow[3] == contigInfo[2] && blastRow[4] == "0" && blastRow[5] == "0") 
      || (blastRow[2] == "100.00" && blastRow[3] == contigInfo[2] && blastRow[4] == "0" && blastRow[5] == "0"))
    {
      blastMap[std::atoi(contigInfo[0].c_str())].push_back(blastRecord(blastRow[0], std::atoi(blastRow[8].c_str()), std::atoi(blastRow[9].c_str())));
    }
  }
}


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

bool comparePathes(VertexPtrVec pathesToFind, pathNode* pNode, bool isForward)
{
  auto currentNode = pNode;
  if ((isForward && pNode->realForward == FORWARD) || (!isForward && pNode->realForward == REVERSE))
  {
    for (int i = 0; i < pathesToFind.size(); ++i)
    {
      if (currentNode->current == pathesToFind[i])
      {
        if (currentNode->forward != 0)
          currentNode = currentNode->forward;
        else
        {
          if (i != pathesToFind.size() - 1)
            break;
        }
      }
      else
      {
        break;
      }
      if (i == pathesToFind.size() - 1)
      {
        return true;
      }
    }
  }
  else
  if ((!isForward && pNode->realReverse == REVERSE) || (isForward && pNode->realReverse == FORWARD))
  {
    for (int i = 0; i < pathesToFind.size(); ++i)
    {
      if (currentNode->current == pathesToFind[i])
      {
        if (currentNode->reverse != 0)
          currentNode = currentNode->reverse;
        else
        {
          if (i != pathesToFind.size() - 1)
            break;
        }
      }
      else
      {
        break;
      }
      if (i == pathesToFind.size() - 1)
      {
        return true;
      }
    }
  }
  return false;
}

std::string tempIndex = "0";
int extendpaths(SGWalkVector& allpaths, bool isForward)
{
  int totalPathSize = 0;
  for (int i = 0; i < allpaths.size(); ++i)
  {

    auto path = allpaths[i];
    auto vertexSeq = path.getVertices();
    auto firstVertex = vertexSeq[0];
    int currPathesSize = 0;
    for (auto d : firstVertex->getPathes())
    {
      if (comparePathes(vertexSeq, d, isForward))
      {
        currPathesSize++;
      }
    }
    totalPathSize += std::max(1, currPathesSize);
  }
  return totalPathSize;
}

enum Contained
{
  FIRSTINSECOND, SECONDINFIRST, NONE
};

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

std::vector<int> differences;

int uncontainedpaths_SGA(Vertex* startID, ScafVector& scaffolds, bool isForward, StringGraph* graph)
{
  int answer = 0;
  SGWalkVector superWalks;
  for (int i = 0; i < scaffolds.size(); ++i)
  {
    SGWalkVector outWalks;
    SGSearch::findWalks(startID, graph->getVertex(scaffolds[i].first), isForward ? ED_SENSE : ED_ANTISENSE, insertSize, 3000, false, outWalks);
    logVerbose("Number of pathes from " + startID->getID() + " to " + scaffolds[i].first + " is " + std::to_string(outWalks.size()));
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
          differences.push_back(outWalks[j].getEndToStartDistance() - scaffolds[i].second);
          for (auto d : startID->getPathes())
          {
            if (comparePathes(outWalks[j].getVertices(), d, isForward))
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

  logVerbose("Total number of uncontained pathes - " + std::to_string(answer));

  for (auto d : CNByNeighboursTemp)
  {
    if (d.second > CNbyNeighbours[d.first])
      CNbyNeighbours[d.first] = d.second;
  }

  int res = extendpaths(allpaths, isForward);
  logVerbose("Number of pathes after extension step - " + std::to_string(res));

  return res;
}

void DeleteAllEdgesExcept(Vertex* from, Vertex* to, SGWalkVector& outwalks)
{
  for (int i = 0; i < outwalks.size(); ++i)
  {
    if (outwalks[i].getLastVertex() != to)
    {
      auto edges = from->getEdges();
      for (int j = 0; j < edges.size(); ++j)
      {
        if (edges[j]->getEnd() == outwalks[i].getLastVertex())
        {
          if (edges[j]->getTwin() != NULL || edges[j]->getTwin() != edges[j])
            edges[j]->getEnd()->deleteEdge(edges[j]->getTwin());
          from->deleteEdge(edges[j]);
        }
      }
    }
  }
}

void improveGraph(StringGraph* graph, SGWalk& walk)
{
  auto vertices = walk.getVertices();
  for (int i = 0; i < vertices.size() - 1; ++i)
  {
    if (vertices[i]->getEstimatedCN() == 1 && vertices[i]->getSeqLen() > 5000)
    {
      SGWalkVector outwalks;
      SGWalkVector outwalks2;

      SGSearch::findAllWalks(vertices[i], ED_SENSE, 1, 2, outwalks);
      SGSearch::findAllWalks(vertices[i], ED_ANTISENSE, 1, 2, outwalks2);

      bool forward = true;
      for (int j = 0; j < outwalks.size(); ++j)
      {
        if (outwalks[j].getLastVertex() == vertices[i + 1])
        {
          forward = false;
          DeleteAllEdgesExcept(vertices[i], vertices[i + 1], outwalks);
          break;
        }
      }

      if (forward)
      {
        DeleteAllEdgesExcept(vertices[i], vertices[i + 1], outwalks2);
      }

    }
  }


  for (int i = vertices.size() - 1; i > 0; --i)
  {
    if (vertices[i]->getEstimatedCN() == 1 && vertices[i]->getSeqLen() > 500)
    {
      SGWalkVector outwalks;
      SGWalkVector outwalks2;

      SGSearch::findAllWalks(vertices[i], ED_SENSE, 1, 2, outwalks);
      SGSearch::findAllWalks(vertices[i], ED_ANTISENSE, 1, 2, outwalks2);

      bool forward = true;
      for (int j = 0; j < outwalks.size(); ++j)
      {
        if (outwalks[j].getLastVertex() == vertices[i - 1])
        {
          forward = false;
          DeleteAllEdgesExcept(vertices[i], vertices[i - 1], outwalks);
          break;
        }
      }

      if (forward)
      {
        DeleteAllEdgesExcept(vertices[i], vertices[i - 1], outwalks2);
      }
    }
  }

}

void playWithpaths(StringGraph* graph, std::vector<StringVector>& scaffolds)
{
  for (int i = 0; i < scaffolds.size(); ++i)
  {
    for (int j = 0; j < scaffolds[i].size() - 1; ++j)
    {
      SGWalkVector outwalks;
      SGWalkVector outwalks2;

      SGSearch::findWalks(graph->getVertex(scaffolds[i][j]), graph->getVertex(scaffolds[i][j + 1]), ED_SENSE, insertSize, 2000, false, outwalks);
      SGSearch::findWalks(graph->getVertex(scaffolds[i][j]), graph->getVertex(scaffolds[i][j + 1]), ED_ANTISENSE, insertSize, 2000, false, outwalks2);

      if (outwalks.size() == 0 && outwalks2.size() == 0)
        continue;
      else
      {
        if (outwalks2.size() == 0)
        {
          improveGraph(graph, outwalks[0]);
          continue;
        }
        if (outwalks.size() == 0)
        {
          improveGraph(graph, outwalks2[0]);
          continue;
        }
//        improveGraph(graph, outwalks[0].getEndToStartDistance() < outwalks2[0].getEndToStartDistance() ? outwalks[0] : outwalks2[0]);

      }
    }
  }
}

void mapPathesToGraph(StringGraph * graph, blastVector& blastVector)
{
  std::sort(blastVector.begin(), blastVector.end(), CompareBlastRec());
  pathNode* start = new pathNode(graph->getVertex(blastVector[0].contigName));
  graph->getVertex(blastVector[0].contigName)->addPath(start);

  for (int i = 0; i < blastVector.size() - 1; ++i)
  {
    if (blastVector[i].end - blastVector[i + 1].start >= 55)
    {
      pathNode * curr = new pathNode(graph->getVertex(blastVector[i + 1].contigName), graph->getVertex(blastVector[i].contigName)->getPathes()[graph->getVertex(blastVector[i].contigName)->getPathes().size() - 1]);
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
      graph->getVertex(blastVector[i].contigName)->setForwardToPath(graph->getVertex(blastVector[i + 1].contigName)->getPathes()[graph->getVertex(blastVector[i + 1].contigName)->getPathes().size() - 1]);
    }
    else
    {
      pathNode * curr = new pathNode(graph->getVertex(blastVector[i + 1].contigName));
      graph->getVertex(blastVector[i + 1].contigName)->addPath(curr);
    }
  }
}

int delete_work(int argc, char* argv[])
{
  std::string filename = argv[2];
  std::string vertexName = argv[3];
  int distance = std::atoi(argv[4]);
  std::string outputGraphName = argv[5];
  StringGraph* graph = SGUtil::loadASQG(filename, 0);
  SGWalkVector outwalks;
  SGWalkVector outwalks2;
  SGSearch::findAllWalks(graph->getVertex(vertexName), ED_SENSE, 1000, 1000, outwalks);
  SGSearch::findAllWalks(graph->getVertex(vertexName), ED_ANTISENSE, 1000, 1000, outwalks2);
  std::unordered_set<std::string> surroundingVertices;
  std::unordered_set<std::string> allVertices;

  for (auto d : graph->getAllVertices())
  {
    allVertices.insert(d->getID());
  }

  for (int i = 0; i < outwalks.size(); ++i)
  {
    for (auto d : outwalks[i].getVertices())
    {
      surroundingVertices.insert(d->getID());
    }
  }

  for (int i = 0; i < outwalks2.size(); ++i)
  {
    for (auto d : outwalks2[i].getVertices())
    {
      surroundingVertices.insert(d->getID());
    }
  }

  for (auto d : allVertices)
  {
    if (surroundingVertices.find(d) == surroundingVertices.end())
    {
      graph->removeConnectedVertex(graph->getVertex(d));
    }
  }
  graph->writeASQG(outputGraphName);
  return 0;
}


int main_work(int argc, char* argv[])
{

  // Arguments
  // primary2-graph.asqg -- SGA graph structure
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
  //  std::cout << blastContigsCommand << std::endl;
 // system(blastContigsCommand.c_str());

  //if (excode)
  //{
  //  std::cerr << "Error making blast db" << std::endl;
  //  return 1;
  //}


  //contruct adjacency lists
  StringGraph* graph = SGUtil::loadASQG(filename, 0);
  blastVector blastOfAllContigs;
  parseBlastFile(blastOfAllContigs, graph);
  mapPathesToGraph(graph, blastOfAllContigs);



  std::ifstream CNEstimation(argv[2]);

  auto allVertices = graph->getAllVertices();

  std::ifstream scafParser(argv[3]);

  parseDEScaffoldFile(scafParser);


  std::ifstream scaffoldFile(argv[4]);

  mapCNtoContigs(graph, CNEstimation);

  //  solveSimpleCases(graph);
//    std::vector<StringVector> scaffolds = parseScaffoldFile(scaffoldFile);
//    playWithpaths(graph, scaffolds);


  for (auto d : graph->getAllVertices())
  {
    std::cout << d->getID() << " is processed" << std::endl;
    std::string idname = d->getID();
    if (idname == "contig-272")
    {
      int z = 0;
    }
    if (d->getSeqLen() > 10000)
    {
      logVerbose(d->getID() + " is skipped. High length.");
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
        isTrulyProcessed = true;

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

  //solveSimpleCases(graph);

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

  //solveSimpleCases(graph);


  // Write the results
  writeCopyNumbers(graph);
  std::sort(differences.begin(), differences.end());
  SGFastaVisitor av("new_contigs.fasta");
  graph->visit(av);

  //  doppelGraph a(graph);
  graph->writeASQG(filename + ".out");

  std::sort(differences.begin(), differences.end());
  for (auto d : differences)
  {
    std::cout << d <<",";
  }

  return 0;
}

int main(int argc, char* argv[])
{
  std::string temp = argv[1];
  if (temp == "delete")
  {
    delete_work(argc, argv);
  }
  else
  {
    main_work(argc, argv);
  }
}

