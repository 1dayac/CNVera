//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// Bidirectional graph 
//
#ifndef BIGRAPH_H
#define BIGRAPH_H

#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include "GraphCommon.h"
#include "Vertex.h"
#include "Edge.h"

//
// Typedefs
//
typedef std::map<VertexID, Vertex*> VertexPtrMap;
typedef std::map<VertexID, Vertex*>::iterator VertexPtrMapIter;
typedef std::map<VertexID, Vertex*>::const_iterator VertexPtrMapConstIter;

//typedef __gnu_cxx::unordered_map<VertexID, Vertex*> VertexPtrMap;
//typedef std::tr1::unordered_map<VertexID, Vertex*> VertexPtrMap;
//typedef SparseHashMap<VertexID, Vertex*, StringHasher> VertexPtrMap;

//typedef VertexPtrMap::iterator VertexPtrMapIter;
//typedef VertexPtrMap::const_iterator VertexPtrMapConstIter;

class Bigraph;
typedef bool(*VertexVisitFunction)(Bigraph*, Vertex*);

typedef EdgePtrVec Path; // alias
typedef std::vector<Path> PathVector;
typedef std::vector<VertexID> VertexIDVec;
typedef std::vector<Vertex*> VertexPtrVec;

class Bigraph
{

    public:
    
        Bigraph();
        ~Bigraph();

        // Add a vertex
        void addVertex(Vertex* pVert);
        
        // Remove a vertex
        // removeIslandVertex removes a vertex that is guarenteed to not have edges
        // removeConnectedVertex removes a (possibly) connected vertex and all the edges to/from it
        void removeIslandVertex(Vertex* pVertex);
        void removeConnectedVertex(Vertex* pVertex);

        // Check if a vertex exists
        bool hasVertex(VertexID id);

        // Get a vertex
        Vertex* getVertex(VertexID id) const;

        // Add an edge
        void addEdge(Vertex* pVertex, Edge* pEdge);

        // Remove an edge
        void removeEdge(const EdgeDesc& ed);

        // Remove all edges marked by color c
        int sweepVertices(GraphColor c);
        int sweepEdges(GraphColor c);

        // Merge vertices
        void mergeVertices(VertexID id1, VertexID id2);

        // Merge vertices that are joined by the specified edge
        void merge(Vertex* pV1, Edge* pEdge);

        // Rename all the vertices in the graph
        void renameVertices(const std::string& prefix = "");

        // Simplify the graph by removing transitive edges
        void simplify();

        // Validate that the graph is sane
        void validate();

        // Flip a given vertex
        void flip(VertexID id);

        // Sort all the vertex adjacency lists
        void sortVertexAdjListsByLen();
        void sortVertexAdjListsByID();

        // Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
        VertexIDVec getNonBranchingVertices() const;

        // Get the linear components of a non-branching graph
        PathVector getLinearComponents();

        // Return all the path of nodes that can be linearally reached from this node
        // The path expands in both directions so the first node in the path is not necessarily the source
        Path constructLinearPath(VertexID id);

        // Reverse a path
        static Path reversePath(const Path& path);

        // Returns an iterator to the first vertex in the graph.
        // If the graph is empty, a null pointer is returned
        Vertex* getFirstVertex() const;

        // Print simple summary statistics to stdout
        void stats() const;
        void printMemSize() const;

        size_t getNumVertices() const;

        // Visit each vertex in the graph and perform the visit function
        bool visit(VertexVisitFunction f);

        //
        VertexPtrVec getAllVertices() const;

        // Append each vertex sequence to the vector of strings
        void getVertexSequences(std::vector<std::string>& outSequences) const;

        // Visit each vertex in the graph and call the visit functor object
        template<typename VF>
        bool visit(VF& vf)
        {
            bool modified = false;
            vf.previsit(this);
            VertexPtrMapConstIter iter = m_vertices.begin(); 
            for(; iter != m_vertices.end(); ++iter)
            {
                modified = vf.visit(this, iter->second) || modified;
            }
            vf.postvisit(this);
            return modified;
        }
        
        // Set the colors for the entire graph
        void setColors(GraphColor c);

        // Check the colors for the entire graph
        bool checkColors(GraphColor c);

        // Get/Set the containment flag
        void setContainmentFlag(bool b);
        bool hasContainment() const;

        //
        void setTransitiveFlag(bool b);
        bool hasTransitive() const;

        //
        void setMinOverlap(int mo);
        int getMinOverlap() const;

        //
        void setErrorRate(double er);
        double getErrorRate() const;

        void setExactMode(bool b);
        bool isExactMode() const;

        // Write the graph to a file
        void writeDot(const std::string& filename, int dotFlags = 0) const;
        void writeASQG(const std::string& filename) const;

        // Returns an allocator for the edges of the graph
        std::set<Edge>* getEdgeAllocator() { return m_pEdgeAllocator; }

        // Returns an allocator for the vertices of the graph
        std::set<Vertex>* getVertexAllocator() { return m_pVertexAllocator; }

        // Return a string for a color code
        static std::string getColorString(GraphColor c);

    private:
        
        // Simplify the graph by compacting edges in the given direction
        void simplify(EdgeDir dir);

        void followLinear(VertexID id, EdgeDir dir, Path& outPath);

        //
        // data
        //
        VertexPtrMap m_vertices;

        // Graph parameters
        bool m_hasContainment;
        bool m_hasTransitive;
        bool m_isExactMode;

        int m_minOverlap;
        double m_errorRate;

        // Memory management
        std::set<Vertex>* m_pVertexAllocator;
        std::set<Edge>* m_pEdgeAllocator;
};

#endif
