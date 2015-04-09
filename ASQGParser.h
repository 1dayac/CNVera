#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

struct asqg_node
{
public:
  unsigned int leftID;
  unsigned int rightID;
  std::vector<short> coverage;
  std::string sequence = "";
  unsigned int numberOfSupportingReads;
  std::unordered_map<int, int> leftNeighbours; //ID of neighbour + overlap
  std::unordered_map<int, int> rightNeighbours; //ID of neighbour + overlap
  asqg_node(unsigned int _leftID, unsigned int _rightID)
    : leftID(_leftID), rightID(_rightID)
  { };
};


struct asqg_graph
{
  std::unordered_map<int, std::shared_ptr<asqg_node>> nodes;
  mag_graph()
  { };
  void readFromFile(std::string filename);
  void outputStats();
  void writeToFile(std::string filename);
};