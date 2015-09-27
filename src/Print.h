/**
 * \file
 * \brief Printout data and information
 */
#ifndef _Print_h
#define _Print_h
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

/**
 * Wrap string by word, by inserting newlines, with maximum length of \p
 * line_width
 *
 * \param[in] string     String to wrap by word
 * \param[in] line_width Maximum length of string before wrapping.
 *
 * \return \p string wrapped by word
 */
char * word_wrap (char* string, int line_width);

/**
 * Indent sections of Help() to aid readability
 *
 * \param[in] text   String to print out
 * \param[in] indent Number of characters to indent \p text
 * \param[in] wrap   Length of string, for word wrapping.
 */
void HelpIndented(const char *text, unsigned int indent, unsigned int wrap);

/**
 * Show help for an option, listing the short, and long-hand notations.
 *
 * \param[in] v1   Short-hand notation of a command line option
 * \param[in] v2   Long-hand notation of a command line option
 * \param[in] text Description of the command line option
 */
void HelpOption(const char *v1, const char *v2, const char *text);

/**
 * Show help
 *
 * \param[in] callname  Name of this program as called from the command line.
 */
void Help(char *callname);

#endif
