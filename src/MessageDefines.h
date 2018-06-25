/**
 * \file
 * \date  28 May 2015
 * \brief Various macros for messages and debugging
 *
 **/
#ifndef _MessageDefines_h
#define _MessageDefines_h
#include <iostream>

// namespace Color {
//     enum Code {
//         FG_RED      = 31,
//         FG_GREEN    = 32,
//         FG_BLUE     = 34,
//         FG_DEFAULT  = 39,
//         BG_RED      = 41,
//         BG_GREEN    = 42,
//         BG_BLUE     = 44,
//         BG_DEFAULT  = 49
//     };
//     std::ostream& operator<<(std::ostream& os, Code code) {
//         return os << "\033[" << static_cast<int>(code) << "m";
//     }
// }

#ifdef DEBUG
/** Macro for writing debug messages, colored for easy viewing.  */
#define DEBUG_MSG(str) do { if (DEBUG) std::cout << "\033[38;5;1mDEBUG: \033[38;5;9m" << str << "\033[39m\n"; } while ( false )
#else
/** Macro for writing debug messages.  */
#define DEBUG_MSG(str) do { } while ( false )
#endif

extern bool THB_VERBOSE;

#ifdef __linux
	/** \name Macros for writing verbose messages to stdout*/
	/**@{*/
	/** Write \p str followed by a newline, \c '\\n'  */
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n"; } while ( false )
	/** Write \p str */
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str; } while ( false )
	/** Write \p str followed by a return, \c '\\r'  */
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
	/**@}*/
	/** \name Macros for writing brief messages to stdout */
	/**@{*/
	/** Write \p str followed by a newline, \c '\\n' */
	#define BRIEF_MSG(str)  do { std::cout << str << "\n"; } while ( false )
	/** Write \p str */
	#define BRIEF_CMSG(str) do { std::cout << str; } while ( false )
	/** Write \p str followed by a return, \c '\\r' */
	#define BRIEF_RMSG(str) do { std::cout << str << "\r" << std::flush; } while ( false )
	/**@}*/
#elif _WIN32
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n" << std::flush; } while ( false )
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str << std::flush; } while ( false )
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
	#define BRIEF_MSG(str)  do { std::cout << str << "\n" << std::flush; } while ( false )
	#define BRIEF_CMSG(str) do { std::cout << str << std::flush; } while ( false )
	#define BRIEF_RMSG(str) do { std::cout << str << "\r" << std::flush; } while ( false )
#endif

#endif // _MessageDefines_h
