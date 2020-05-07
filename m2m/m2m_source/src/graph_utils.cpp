/******************************************************************************
* A class for analyzing the map topology (graphs).
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#include <list>
#include <exception>

#include "graph_utils.h"

namespace ConverterLib {

int RemoveUncommonUnits(Graph &A, const Graph &B)
{
  int num_edges_removed = 0;

  // Loop over all edges of graph A.
  UnitContacts::iterator A_it = A.begin();
  while (A_it != A.end())
  {
    bool unit1_found = false;
    bool unit2_found = false;

    for (UnitContacts::const_iterator B_it = B.begin(); B_it != B.end(); ++B_it)
    {
      if (B_it->unit1->name == A_it->unit1->name || B_it->unit2->name == A_it->unit1->name) unit1_found = true;
      if (B_it->unit1->name == A_it->unit2->name || B_it->unit2->name == A_it->unit2->name) unit2_found = true;

      if (unit1_found && unit2_found) break;
    }

    if (unit1_found && unit2_found)
    {
      A_it++;
    }
    else
    {
      // Remove this edge from the graph.
      A_it = A.erase(A_it);
      num_edges_removed++;
    }
  }
  return num_edges_removed;
}

//=======================================================================================================
// indexes[0] - normal index,
// indexes[1] - weighted index (by faultiness).
void CalculateJaccardIndex(const Graph &A, const Graph &B, std::vector<double> &indexes)
{
  indexes.resize(2);

  if (A.size() == 0 && B.size() == 0)
  {
    // By definition.
    indexes[0] = 1.0;
    indexes[1] = 1.0;
  }
  else
  {
    size_t intersection_number = 0;
    double intersection_number_weighted = 0.0;

    // Loops that are comparing all pairs of edges between graphs A and B.
    for (UnitContacts::const_iterator A_it = A.begin(); A_it != A.end(); ++A_it)
      for (UnitContacts::const_iterator B_it = B.begin(); B_it != B.end(); ++B_it)
    {
      if ((A_it->unit1->name == B_it->unit1->name && A_it->unit2->name == B_it->unit2->name)
          || (A_it->unit1->name == B_it->unit2->name && A_it->unit2->name == B_it->unit1->name))
      { // Same graph edge found.
        intersection_number++;

        // Determine weighted intersection.
        double weight1 = A_it->faults_length / A_it->total_length;
        double weight2 = B_it->faults_length / B_it->total_length;

        intersection_number_weighted += 1.0 - std::fabs(weight1 - weight2);
      }
    }

    int denom = (int)(A.size() + B.size() - intersection_number);
    double denom_weighted = double(A.size() + B.size()) - intersection_number_weighted;

    if (denom == 0 || denom_weighted == 0) {
      throw std::runtime_error("Error: wrong denominator in CalculateJaccardIndex()!");
    }

    indexes[0] = (double(intersection_number) / double(denom));
    indexes[1] = (intersection_number_weighted / denom_weighted);
  }
}

//=======================================================================================
double CalculateGraphDensity(const Graph &graph)
{
  // A list of graph nodes.
  std::list<std::string> graph_nodes;
  // Building a list of nodes.
  for (Graph::const_iterator it = graph.begin(); it != graph.end(); ++it)
  {
    graph_nodes.push_back(it->unit1->name);
    graph_nodes.push_back(it->unit2->name);
  }
  graph_nodes.sort();
  graph_nodes.unique();

  size_t num_nodes = graph_nodes.size();
  size_t num_edges = graph.size();

  double density;

  if (num_nodes > 0)
  {
    density = double(2 * num_edges) / double(num_nodes * (num_nodes - 1));
  }
  else
  {
    // Define it this way.
    density = 1.0;
  }

  return density;
}

}

