#include "ASQGParser.h"
#include <fstream>  
#include <sstream>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <map>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

void asqg_graph::readFromFile(std::string filename)
{
  std::ifstream in(filename, std::ios::in);
  if (!in.is_open())
  {
    std::cerr << "Error! File couldn't be opened!" << std::endl;
    return;
  }
  while (!in.eof())
  {
    std::string str = "";
    std::getline(in, str);
    if (str.substr(0, 2) == "HT")
    {
      parseHeader(str);
    }
    else
    if (str.substr(0, 2) == "VT")
    {
      parseVertex(str)
    }
    else
    if (str.substr(0, 2) == "ED")
    {
      parseEdge(str)
    }
    else
    {
      continue;
    }
  }
}

void asqg_graph::outputStats()
{
  //TODO
  int consistentCoverages = 0;
  int unconsistentCoverages = 0;
  std::map<int, int> nodeDegree;
  for (auto it : nodes)
  {
    nodeDegree[it.second->leftNeighbours.size() + it.second->rightNeighbours.size()]++;
    if (it.second->leftNeighbours.size() + it.second->rightNeighbours.size() >= 5)
    {
      std::cout << it.second->leftID << ":" << it.second->rightID << std::endl;
    }
    if (it.second->coverage.size() > 220)
    {
      auto result = std::minmax_element(it.second->coverage.begin() + 110, it.second->coverage.end() - 110);
      if (*result.second - *result.first > 20)
      {
        unconsistentCoverages++;
      }
      else
      {
        consistentCoverages++;
      }
    }
  }
  std::cout << "Number of nodes with unconsistent coverage - " << unconsistentCoverages << std::endl;
  std::cout << "Number of nodes with consistent coverage - " << consistentCoverages << std::endl;
  std::cout << "Summary about node degrees:" << std::endl;
  for (auto it : nodeDegree)
  {
    it.second /= 2;
    std::cout << it.first << " - " << it.second << std::endl;
  }
}

void asqg_graph::writeToFile(std::string filename)
{
  std::ofstream out(filename, std::ios::out);
  std::unordered_set<int> used_nodes;
  for (auto it : nodes)
  {
    if (used_nodes.find(it.first) != used_nodes.end())
    {
      continue;
    }
    used_nodes.insert(it.second->leftID);
    used_nodes.insert(it.second->rightID);
    out << '@' << it.second->leftID << ':' << it.second->rightID << '\t';
    out << it.second->numberOfSupportingReads << '\t';
    for (auto it2 : it.second->leftNeighbours)
    {
      out << it2.first << ',' << it2.second << ';';
    }
    out << '\t';
    for (auto it2 : it.second->rightNeighbours)
    {
      out << it2.first << ',' << it2.second << ';';
    }
    out << '\n' << it.second->sequence << "\n+\n";
    for (auto it2 : it.second->coverage)
      out << it2;
    out << '\n';
  }
}