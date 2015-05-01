#ifndef _NeighborPrint_h
#define _NeighborPrint_h
#include <vector>
#include <math.h>
#include "Histograms.h"
#include "SimpleMath.h"

// void Print_XYZ( std::vector<double> x,
//                 std::vector<double> y,
//                 std::vector<double> z);

void Print_AllFrames(std::ostream *out,
                     std::vector<struct Histograms_s> *frame);
void Print_CombineFrames(std::ostream *out,
                         std::vector<struct Histograms_s> *frame);
void Print_CombineNeighbors(std::ostream *out,
                            std::vector<struct Histograms_s> *frame);
#endif // NeighborPrint
// vim:tw=76:ts=4:sw=4
