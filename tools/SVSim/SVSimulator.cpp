// SVSimulator.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fasta.h"

void usage()
{
  std::cerr << "This script takes as first parameter file with multiple lines, each line is in following format:" << std::endl;
  std::cerr << "X1 X2 X3" << std::endl;
  std::cerr << "X1 = 1 if copy and X1 = 2 if delete" << std::endl;
  std::cerr << "X2 - length of copied or deleted fragment" << std::endl;
  std::cerr << "X3 - number of copies, if X1 = 1" << std::endl;
  std::cerr << "X3 - number of copies, if X1 = 1" << std::endl;
  std::cerr << "This script takes as second parameter file with reference in fasta format" << std::endl;
}


int main(int argc, char* argv[])
{

  if (argc != 3)
  {
    usage();
    return 1;
  }
  srand(time(NULL));

  std::vector<FastaRecord> record;
  FastaReader* f = new FastaReader(argv[2]);
  f->GetSequences(record);

  std::ifstream inSV(argv[1]);
  
  std::ofstream logs("logs.txt");
  int type = 0;
  while (inSV >> type)
  {
    int len = 0;
    int num = 0;
    if (type == 1)
    {
      inSV >> len >> num;
      int seq = rand() % record.size();
      int pos = rand() % record[seq].sequence_.length();
      std::string stringToInsert = record[seq].sequence_.substr(pos, len);
      logs << "Inserted " << len << " from " << record[seq].id_ << " to:" << std::endl;

      for (int i = 0; i < num; ++i)
      {

        int seq2 = rand() % record.size();
        int pos2 = rand() % record[seq].sequence_.length();
        record[seq].sequence_ = record[seq].sequence_.substr(0, pos2) + stringToInsert + record[seq].sequence_.substr(pos2, record[seq].sequence_.length() - pos2);

      }
    }
    else
    if (type == 2)
    {
      inSV >> len;
      int seq = rand() % record.size();
      int pos = rand() % record[seq].sequence_.length();
      record[seq].sequence_.erase(pos, len);
      logs << "Deleted " << len << " from " << record[seq].id_ << " at " << pos << std::endl;
    }
    else
    {
      continue;
    }
  }
  std::string outFilename = argv[2];

  FastaWriter out(outFilename.substr(0, outFilename.find_last_of('.')) + ".sim" + outFilename.substr(outFilename.find_last_of('.'), outFilename.length() - outFilename.find_last_of('.') + 1));
  out.writeSequences(record);
	return 0;
}

