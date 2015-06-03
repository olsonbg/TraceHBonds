#ifndef _MessageDefines_h
#define _MessageDefines_h
#include <ostream>

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
#define DEBUG_MSG(str) do { if (DEBUG) std::cout << "\033[38;5;1mDEBUG: \033[38;5;9m" << str << "\033[39m\n"; } while ( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#ifdef __linux
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n"; } while ( false )
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str; } while ( false )
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
	#define BRIEF_MSG(str)  do { std::cout << str << "\n"; } while ( false )
	#define BRIEF_CMSG(str) do { std::cout << str; } while ( false )
	#define BRIEF_RMSG(str) do { std::cout << str << "\r" << std::flush; } while ( false )
#elif _WIN32
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n" << std::flush; } while ( false )
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str << std::flush; } while ( false )
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
	#define BRIEF_MSG(str)  do { std::cout << str << "\n" << std::flush; } while ( false )
	#define BRIEF_CMSG(str) do { std::cout << str << std::flush; } while ( false )
	#define BRIEF_RMSG(str) do { std::cout << str << "\r" << std::flush; } while ( false )
#endif

#endif // _MessageDefines_h
