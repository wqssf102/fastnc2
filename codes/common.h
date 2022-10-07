#ifndef __COMMON_H__
#define __COMMON_H__


#include "local_config.h"


#include "fastnc.h"

// Convert character to integer (for commandline argument parsing)
int int_from_optarg(const char *optarg);

// Convert character to float (for commandline argument parsing)
float float_from_optarg(const char *optarg);

#endif