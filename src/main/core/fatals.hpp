#pragma once

#include "definesStandard.hpp"
#include "logger.hpp"

/************************************************************************************************************************ 
 *  @brief Checks a boolean expression. If false, a fatal-error is passed to the prompt. Cannot be disabled.
 * 
 *  @param expr    a boolean expression to check (e.g. 1 == 0)
 *  @param message a message to throw before the code crashs. Message uses format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define FATAL_ASSERT(expr, message)                                                                                    \
{                                                                                                                      \
    if(expr){}                                                                                                         \
    else {                                                                                                             \
        FATAL_MSG("Assertion Failure: %s, message: %s, in file: %s, line: %d\n", #expr, message, __FILE__, __LINE__);  \
        exit(EXIT_FAILURE_ASSERTION);                                                                                                       \
    }                                                                                                                  \
}                                                                                                                      \


