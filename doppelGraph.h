#include "SGUtil.h"
#include <unordered_map>
//monotization of bidirected graph
//and Dinic algorithm on it
struct doppelGraph
{
public:
  doppelGraph();
  doppelGraph(StringGraph* graph);
private:
  std::vector<std::vector<int>> incidenceMatrix;
  std::vector<std::vector<int>> flowMatrix;
  std::unordered_map<int, std::string> vertexLabels;
  std::unordered_map<std::string, int> reverseVertexLabels;

  
};