#include "common.h"
//

int int_from_optarg(const char *optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(optarg);
    std::string string_int = optstring.substr(0, 8);
    for (std::string::iterator it = string_int.begin(); it != string_int.end(); ++it) {
        if (!isdigit(*it)) {
            fprintf(stderr, "This doesn't look like a usable integer: %s\n", optarg);
            exit(1);
        }
    }
    return std::atoi(string_int.c_str());
}


float float_from_optarg(const char *optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(optarg);
    std::string string_float = optstring.substr(0, 8);
    for (std::string::iterator it = string_float.begin(); it != string_float.end(); ++it) {
        if (!isdigit(*it) && (*it) != '.') {
            fprintf(stderr, "This doesn't look like a usable float: %s\n", optarg);
            exit(1);
        }
    }
    return std::atof(string_float.c_str());
}