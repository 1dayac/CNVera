#include "doppelGraph.h"

doppelGraph::doppelGraph()
{ }

doppelGraph::doppelGraph(StringGraph* strGraph)
{
  VertexPtrVec vertexVec = strGraph->getAllVertices();

  for (int i = 0; i < vertexVec.size(); ++i)
  {
    vertexLabels[i] = vertexVec[i]->getID();
    vertexLabels[i + vertexVec.size()] = vertexVec[i]->getID() + "_2";
    reverseVertexLabels[vertexVec[i]->getID()] = i;
    reverseVertexLabels[vertexVec[i]->getID() + "_2"] = i + vertexVec.size();
  }
  incidenceMatrix.resize(vertexVec.size()*2);
  for (int i = 0; i < incidenceMatrix.size(); ++i)
  {
    incidenceMatrix[i].resize(vertexVec.size() * 2);
  }

  for (int i = 0; i < vertexVec.size(); ++i)
  {
    EdgePtrVec edges = vertexVec[i]->getEdges();
    for (int j = 0; j < edges.size(); ++j)
    {
      EdgeDir dir = edges[j]->getDir();
      EdgeComp comp = edges[j]->getComp();
      if (edges[j]->getDesc().pVertex->getID() == vertexVec[i]->getID())
      {
        //loop
        if (edges[j]->getDir() == ED_SENSE)
        {
          incidenceMatrix[i][i + vertexVec.size()] = 2;
        }
        else
        {
          incidenceMatrix[i + vertexVec.size()][i] = 2;
        }
      }
      else
      if ((edges[j]->getDir() == ED_SENSE && edges[j]->getComp() == EC_REVERSE) || (edges[j]->getDir() == ED_ANTISENSE && edges[j]->getComp() == EC_SAME))
      {
        //different signs
        if (edges[j]->getDir() == ED_SENSE)
        {
          incidenceMatrix[i][reverseVertexLabels[edges[j]->getDesc().pVertex->getID()]] = 1;
          incidenceMatrix[reverseVertexLabels[edges[j]->getDesc().pVertex->getID()] + vertexVec.size()][i + vertexVec.size()] = 1;
        }
        else
        {
          incidenceMatrix[reverseVertexLabels[edges[j]->getDesc().pVertex->getID()]][i] = 1;
          incidenceMatrix[i + vertexVec.size()][reverseVertexLabels[edges[j]->getDesc().pVertex->getID()] + vertexVec.size()] = 1;
        }
      }
      else
      {
        //one sign
        if (edges[j]->getDir() == ED_SENSE)
        {
          incidenceMatrix[i][reverseVertexLabels[edges[j]->getDesc().pVertex->getID()] + vertexVec.size()] = 1;
          incidenceMatrix[reverseVertexLabels[edges[j]->getDesc().pVertex->getID()]][i + vertexVec.size()] = 1;
        }
        else
        {
          incidenceMatrix[reverseVertexLabels[edges[j]->getDesc().pVertex->getID()] + vertexVec.size()][i] = 1;
          incidenceMatrix[i + vertexVec.size()][reverseVertexLabels[edges[j]->getDesc().pVertex->getID()]] = 1;
        }
      }
    }
  }
  
}
