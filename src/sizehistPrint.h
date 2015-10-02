/**
 * \file
 * \brief Printout data and information
 */
#ifndef _sizehistPrint_h
#define sizehistPrint_h
#include <iomanip>
#include "TraceHBondsConfig.h"
#include "TraceHBonds.h"

/**
 * Round to the nearest integer, 0.5 gets rounded up.
 *
 * \param[in] v double
 *
 * \return \p v rounded to the nearest integer
 */
inline int Round( double v);

/**
 *
 * Print out the histogram for hydrogen bond chain lengths
 *
 * \param[in] out        Stream to send output to
 * \param[in] Histogram  Histogram data.
 * \param[in] Stop       Last bin of histogram to show
 * \param[in] Start      First bin of historgram to show
 * \param[in] Step       Show every \p step bin
 * \param[in] BarLength  Maximum length of histogram bar, in characters.
 * \param[in] NumBins    Number of bins to show, all from \c Start to \c Stop
 *                       if \c zero
 * \param[in] CC         String to denote comments in output
 *
 */
void PrintHistogramChain( std::ostream *out,
                          std::vector<unsigned int>Histogram,
                          unsigned int Stop, unsigned int Start,
                          unsigned int Step, int BarLength,
                          unsigned int NumBins,
                          std::string CC );

/**
 *
 * Print out the histogram for molecules in hydrogen bond chains.
 *
 * \param[in] out        Stream to send output to
 * \param[in] Histogram  Histogram data.
 * \param[in] Max        Last bin of histogram to show
 * \param[in] Min        First bin of historgram to show
 * \param[in] Step       Show every \p step bin
 * \param[in] BarLength  Maximum length of histogram bar, in characters.
 * \param[in] NumBins    Number of bins to show, all from \c Start to \c Stop
 *                       if \c zero
 * \param[in] CC         String to denote comments in output
 *
 */
void PrintHistogramMolecules( std::ostream *out,
                              std::vector<unsigned int>Histogram,
                              unsigned int Max, unsigned int Min,
                              unsigned int Step, int BarLength,
                              unsigned int NumBins,
                              std::string CC );
#endif
