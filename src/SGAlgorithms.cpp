//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#include "SGAlgorithms.h"
#include "SGUtil.h"
#include <iterator>

// add edges to the graph for the given overlap
Edge* SGAlgorithms::createEdgesFromOverlap(StringGraph* pGraph, const Overlap& o, bool allowContained, size_t maxEdges)
{
    // Initialize data and perform checks
    Vertex* pVerts[2];
    EdgeComp comp = (o.match.isRC()) ? EC_REVERSE : EC_SAME;

    bool isContainment = o.match.isContainment();
//    assert(allowContained || !isContainment);
    (void)allowContained;
    for(size_t idx = 0; idx < 2; ++idx)
    {
        pVerts[idx] = pGraph->getVertex(o.id[idx]);

        // If one of the vertices is not in the graph, skip this edge
        // This can occur if one of the verts is a strict substring of some other vertex so it will
        // never be added to the graph
        if(pVerts[idx] == NULL)
            return NULL;
    }

    // Check if this is a substring containment, if so mark the contained read
    // but do not create edges
    for(size_t idx = 0; idx < 2; ++idx)
    {
        if(!o.match.coord[idx].isExtreme())
        {
            size_t containedIdx = 1 - idx;
            assert(o.match.coord[containedIdx].isExtreme());
            pVerts[containedIdx]->setColor(GC_RED);
            pGraph->setContainmentFlag(true);
            return NULL;
        }
    }

    // If either vertex has the maximum number of edges,
    // do not add any more. This is to protect against ultra-dense
    // regions of the graph inflating memory usage. The nodes that reach
    // this limit, and nodes connected to them are marked as super repeats.
    // After loading the graph, all edges to super repeats are cut to prevent
    // misassembly.
    size_t num_edges_0 = pVerts[0]->countEdges();
    size_t num_edges_1 = pVerts[1]->countEdges();
    if(num_edges_0 > maxEdges || num_edges_1 > maxEdges)
    {
        WARN_ONCE("Edge limit reached for vertex when loading graph");
        pVerts[0]->setSuperRepeat(true);
        pVerts[1]->setSuperRepeat(true);
        return NULL;
    }

    if(!isContainment)
    {
        Edge* pEdges[2];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            EdgeDir dir = o.match.coord[idx].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE;
            const SeqCoord& coord = o.match.coord[idx];
            pEdges[idx] = new Edge(pVerts[1 - idx], dir, comp, coord);
        }

        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);
        
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[1], pEdges[1]);
        return pEdges[0];
    }
    else
    {
        // Contained edges don't have a direction, they can be travelled from
        // one vertex to the other in either direction. Hence, we 
        // add two edges per vertex. Later during the contain removal
        // algorithm this is important to determine transitivity
        Edge* pEdges[4];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            const SeqCoord& coord = o.match.coord[idx];
            pEdges[idx] = new(pGraph->getEdgeAllocator()) Edge(pVerts[1 - idx], ED_SENSE, comp, coord);
            pEdges[idx + 2] = new(pGraph->getEdgeAllocator()) Edge(pVerts[1 - idx], ED_ANTISENSE, comp, coord);
        }
        
        // Twin the edges and add them to the graph
        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);

        pEdges[2]->setTwin(pEdges[3]);
        pEdges[3]->setTwin(pEdges[2]);
    
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[0], pEdges[2]);

        pGraph->addEdge(pVerts[1], pEdges[1]);
        pGraph->addEdge(pVerts[1], pEdges[3]);
        
        // Set containment flags
        updateContainFlags(pGraph, pVerts[0], pEdges[0]->getDesc(), o);
        return pEdges[0];
    }
}

//
void SGAlgorithms::updateContainFlags(StringGraph* pGraph, Vertex* pVertex, const EdgeDesc& ed, const Overlap& ovr)
{
    assert(ovr.isContainment());
    // Determine which of the two vertices is contained
    Vertex* pOther = ed.pVertex;
    if(ovr.getContainedIdx() == 0)
        pVertex->setContained(true);
    else
        pOther->setContained(true);
    pGraph->setContainmentFlag(true);
}

