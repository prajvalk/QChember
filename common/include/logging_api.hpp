#ifndef NEWSCF_LOGGING_API
#define NEWSCF_LOGGING_API

#define NEWSCF_ENABLE_LOGGING

#ifdef NEWSCF_ENABLE_LOGGING
#include <string>
#include <iostream>

#ifndef NEWSCF_LOG_LEVEL
#define NEWSCF_LOG_LEVEL WARN
#endif

#define LOG(LEVEL, MSG) \
    newscf_log(LEVEL, MSG, __FILE__, __LINE__)

#define LOG_ASIS(LEVEL, MSG, F, L) \
    newscf_log(LEVEL, MSG, F, L)

enum NEWSCF_LOGLEVELS {
    DEV_DUMP,
    DEV_INFO,
    DEV_WARN,
    DEV_ERROR,
    INFO,
    WARN,
    ERROR
};

inline void newscf_log (NEWSCF_LOGLEVELS level, std::string msg, std::string file, int line) {
    if (level < NEWSCF_LOG_LEVEL) return;

    std::string level_str = "";

    if (level == DEV_DUMP)              level_str = "DEV_DUMP ";
    else if (level == DEV_INFO)         level_str = "DEV_INFO ";
    else if (level == DEV_WARN)         level_str = "DEV_WARN ";
    else if (level == DEV_ERROR)        level_str = "DEV_ERROR";
    else if (level == INFO)             level_str = "INFO     ";
    else if (level == WARN)             level_str = "WARN     ";
    else                                level_str = "ERROR    ";

    std::cout << "==== "+level_str+" ==== ("+file+":"+std::to_string(line)+"): "+msg << std::endl;
}

#else
#define LOG(LEVEL, MSG)
#endif

#endif // NEWSCF_LOGGING_API