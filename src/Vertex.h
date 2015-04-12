//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#ifndef VERTEX_H
#define VERTEX_H

// Includes
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"
#include "QualityVector.h"
#include "EdgeDesc.h"
#include "MultiOverlap.h"

// Forward declare
class Edge;

// Default edge sorting function, by ID
struct EdgeIDComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

// Edge sorting function, by length
struct EdgeLenComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

enum realDirection
{
  NONEDir, FORWARD, REVERSE
};

struct pathNode
{
public:
  pathNode* forward;
  realDirection realForward;
  pathNode* reverse;
  realDirection realReverse;
  Vertex * current;
  pathNode(Vertex * _current, pathNode * _reverse = 0, pathNode* _forward = 0, realDirection _realForward = NONEDir, realDirection _realReverse = NONEDir)
    : current(_current), forward(_forward), reverse(_reverse), realForward(_realForward), realReverse(_realReverse)
  {  }

};


// Typedefs
typedef std::map<EdgeDesc, Edge*> EdgePtrMap;
typedef std::vector<Edge*> EdgePtrVec;
typedef std::set<EdgeDesc> EdgeDescSet;
typedef std::list<Edge*> EdgePtrList;
typedef EdgePtrMap::iterator EdgePtrMapIter;
typedef EdgePtrMap::const_iterator EdgePtrMapConstIter;
typedef EdgePtrVec::iterator EdgePtrVecIter;
typedef EdgePtrVec::const_iterator EdgePtrVecConstIter;
typedef EdgePtrList::iterator EdgePtrListIter;
typedef EdgePtrList::const_iterator EdgePtrListConstIter;
typedef std::vector<pathNode*> pathStack;

class Vertex
{
    public:
    
        Vertex(VertexID id, const std::string& s) : m_id(id), 
                                                    m_seq(s), 
                                                    m_color(GC_WHITE),
                                                    m_coverage(1),
                                                    m_isContained(false),
                                                    m_estimatedCN(0),
                                                    m_estimatedWithNeighboursCN(0),
                                                    m_isTrusted(false),
                                                    m_estimatedWithPairedReads(false)
                                                    {}
        ~Vertex();

        // High-level modification functions
        
        // Merge another vertex into this vertex, as specified by pEdge
        void merge(Edge* pEdge);

        // sort the edges by the ID of the vertex they point to
        void sortAdjListByID();

        // sort the edges by the length of the label of the edge
        void sortAdjListByLen();

        // Ensure that all the edges are unique
        bool markDuplicateEdges(GraphColor dupColor); 

        // Get a multioverlap object representing the overlaps for this vertex
        MultiOverlap getMultiOverlap() const;

        // Edge list operations
        void addEdge(Edge* ep);
        void removeEdge(Edge* pEdge);
        void removeEdge(const EdgeDesc& ed);
        void deleteEdge(Edge* pEdge);
        void deleteEdges();

        int sweepEdges(GraphColor c);
        bool hasEdge(Edge* pEdge) const;
        bool hasEdge(const EdgeDesc& ed) const;
        bool hasEdgeTo(const Vertex* pY) const;

        Edge* getEdge(const EdgeDesc& ed);
        EdgePtrVec findEdgesTo(VertexID id);
        EdgePtrVec getEdges(EdgeDir dir) const;
        EdgePtrVec getEdges() const;
        EdgePtrVecIter findEdge(const EdgeDesc& ed);
        EdgePtrVecConstIter findEdge(const EdgeDesc& ed) const;
        Edge* getLongestOverlapEdge(EdgeDir dir) const;

        size_t countEdges() const;
        size_t countEdges(EdgeDir dir);

        // Calculate the difference in overlap lengths between
        // the longest and second longest edge
        int getOverlapLengthDiff(EdgeDir dir) const;

        // Ensure the vertex data is sane
        void validate() const;
        
        // setters
        void setID(VertexID id) { m_id = id; }
        void setEdgeColors(GraphColor c);
        void setSeq(const std::string& s) { m_seq = s; }
        void setColor(GraphColor c) { m_color = c; }
        void setContained(bool c) { m_isContained = c; }
        void setSuperRepeat(bool b) { m_isSuperRepeat = b; }
        void setEstimatedCN(size_t n) { m_estimatedCN = n; }
        void setTrustedCN(bool b) { m_isTrusted = b; }
        void setEstimatedWithNeighboursCN(size_t n)  { m_estimatedWithNeighboursCN = n;}
        void setEstimatedWithPairedReads(bool b) { m_estimatedWithPairedReads = b; }
        // getters
        VertexID getID() const { return m_id; }
        GraphColor getColor() const { return m_color; }
        const std::string& getSeq() const { return m_seq; }
        std::string getStr() const { return m_seq; }
        size_t getSeqLen() const { return m_seq.length(); }
        bool isContained() const { return m_isContained; }
        bool isSuperRepeat() const { return m_isSuperRepeat; }
        uint16_t getCoverage() const { return m_coverage; }
        size_t getEstimatedCN() const { return m_estimatedCN; }
        size_t getEstimatedWithNeighboursCN() const { return m_estimatedWithNeighboursCN; }
        bool getTrustedCN() const { return m_isTrusted; }
        size_t getUntrustyMetric() { return m_estimatedCN == 0 ? -1 : getSeqLen() / m_estimatedCN; }
        bool getEstimatedWithPairedReads() const { return m_estimatedWithPairedReads; }
        pathStack getPaths() const { return pStack; }
        void setForwardToPath(pathNode* p) { pStack[pStack.size() - 1]->forward = p; }
        void addPath(pathNode* pNode) { pStack.push_back(pNode); };
        void operator delete(void* /*target*/, size_t /*size*/)
        {
            // delete does nothing since all allocations go through the memory pool
            // belonging to the graph. The memory allocated for the vertex will be
            // cleaned up when the graph is destroyed.
        }

        // Output edges in graphviz format
        void writeEdges(std::ostream& out, int dotFlags) const;

    private:

        // Global new is disallowed, all allocations must go through the pool

        // Ensure all the edges in DIR are unique
        bool markDuplicateEdges(EdgeDir dir, GraphColor dupColor);

        VertexID m_id;
        EdgePtrVec m_edges;
        std::string m_seq;
        GraphColor m_color;
        pathStack pStack;
        // Counter of the number of vertices that have been merged into this one
        uint16_t m_coverage; 

        bool m_isContained;
        bool m_isSuperRepeat;

        //Is node trusted (as function of length)?
        bool m_isTrusted;
        bool m_estimatedWithPairedReads;

        size_t m_estimatedCN;
        size_t m_estimatedWithNeighboursCN;
        size_t m_unrustyMetrics;

};

#endif
