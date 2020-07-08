/******************************************************************************
* Functions to work with graphs.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef graph_utils_h
#define graph_utils_h

#include <vector>

#include "converter_types.h"

namespace ConverterLib {

typedef UnitContacts Graph;

//! Calculate Jaccard indexes (for comparing two graphs), normal and weighted.
void CalculateJaccardIndex(const Graph &A, const Graph &B, std::vector<double> &indexes);

//! Calculate the graph density.
double CalculateGraphDensity(const Graph &graph);

//! Remove from graph A nodes (and edges containing those nodes) that are not present in graph B.
//! (E.g. for graphs A=n1-n2-n3, B=n2-n3-n4, the edge n1-n2 in A will be removed, because B has no node n1.)
int RemoveUncommonUnits(Graph &A, const Graph &B);

}

#endif
