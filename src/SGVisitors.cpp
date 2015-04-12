//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGVisitors.h"
//#include "CompleteOverlapSet.h"
//#include "SGSearch.h"
//#include "stdaln.h"

//
// SGFastaVisitor - output the vertices in the graph in 
// fasta format
//
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    m_fileHandle << ">" << pVertex->getID() << " " <<  pVertex->getSeq().length() 
                 << " " << 0 << "\n";
    m_fileHandle << pVertex->getSeq() << "\n";
    return false;
}

