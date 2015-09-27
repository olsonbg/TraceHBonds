/**
 * \file
 * \brief Neighbor histograms
 *
 **/
#ifndef _NeighborPrint_h
#define _NeighborPrint_h
#include <vector>
#include <math.h>
#include "Histograms.h"
#include "SimpleMath.h"

/**
 *
 * Neighbor distances for each frame and all chain lengths
 *
 * \param[in] out    Stream to send results to
 * \param[in] frame  Histograms of neighbor distances of non hydrogen atoms,
 *                   from getNeighbors()
 *
 **/
void Print_AllFrames(std::ostream *out,
                     std::vector<struct Histograms_s> *frame);

/**
 *
 * Combine all frames for an overall average of Neighbor distances
 *
 * \param[in] out    Stream to send results to
 * \param[in] frame  Histograms of neighbor distances of non hydrogen atoms,
 *                   from getNeighbors()
 *
 **/
void Print_CombineFrames(std::ostream *out,
                         std::vector<struct Histograms_s> *frame);
/**
 *
 * Combine all frames and all Nth nearest chain lengths of Neighbor
 * distances
 *
 * \param[in] out    Stream to send results to
 * \param[in] frame  Histograms of neighbor distances of non hydrogen atoms,
 *                   from getNeighbors()
 *
 **/
void Print_CombineNeighbors(std::ostream *out,
                            std::vector<struct Histograms_s> *frame);
#endif // NeighborPrint
// vim:tw=76:ts=4:sw=4
