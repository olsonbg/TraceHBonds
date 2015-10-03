/**
 * \file
 * \brief Printout data and information
 */
#ifndef _help_h
#define _help_h
#include "TraceHBondsConfig.h"
#include "TraceHBonds.h"

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
