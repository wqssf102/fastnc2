#include "fastnc_opts.h"


void print_help() {
    fprintf(stderr, "Program: FastNC (c++ calculate)\n");
    fprintf(stderr, "Contact: Qiusheng WU (565715597@qq.com)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  fastncn [options] --adj_table <file> --outfile <file>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c <file>, --adj_table <file>\n");
    fprintf(stderr, "                The result from get.adjacency function of igraph package\n");
    fprintf(stderr, "  -o <file>, -outfile <file>\n");
    fprintf(stderr, "                Result output table\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -t <float>, -threshold <float>\n");
    fprintf(stderr, "                The threshold for deletion of node (default: 0.8)\n");
	fprintf(stderr, "  -n <int>, -number <int>\n");
    fprintf(stderr, "                Number of iterations (default: 1000)\n");
	fprintf(stderr, "  -j <int>, -job <int>\n");
    fprintf(stderr, "                Number of jobs (default: 4)\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -h        --help\n");
    fprintf(stderr, "                Display this help and exit\n");
}


//
FastNCOptions get_commandline_arguments(int argc, char **argv){
    FastNCOptions fastnc_options;
    struct option long_options[] =
        {
            {"adj_table", required_argument, NULL, 'c'},
            {"outfile", required_argument, NULL, 'o'},
			{"threshold", required_argument, NULL, 't'},
			{"number", required_argument, NULL, 'n'},
			{"job", required_argument, NULL, 'j'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Parse commandline arguments
    while (1) {
        // Parser variables
        int option_index = 0;
        int c;

        // Parser
        c = getopt_long(argc, argv, "hc:o:t:n:j:", long_options, &option_index);

        // If no more arguments to parse, break
        if (c == -1) {
            break;
        }

        // Process current arguments
        switch(c) {
            case 'c':
                fastnc_options.adj_table_filename = optarg;
                break;
            case 'o':
                fastnc_options.outfile_filename = optarg;
                break;
			case 't':
                fastnc_options.threshold_name = float_from_optarg(optarg);
                break;
			case 'n':
                fastnc_options.number_name = int_from_optarg(optarg);
                break;
			case 'j':
                fastnc_options.jobs_name = int_from_optarg(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            default:
                exit(1);
        }
    }


    // Check if have an attempt at arguments
    if (argc < 2) {
        print_help();
        fprintf(stderr,"\n%s: error: option -c/--adj_table and -o/--outfile are required\n", argv[0]);
        exit(1);
    }

    // Abort execution if given unknown arguments
    if (optind < argc){
        print_help();
        fprintf(stderr, "\n%s: invalid argument: %s\n", argv[0], argv[optind++]);
    }


    // Make sure we have filenames
    if (fastnc_options.adj_table_filename.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -c/--adj_table is required\n", argv[0]);
        exit(1);
    }
    if (fastnc_options.outfile_filename.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -o/--outfile is required\n", argv[0]);
        exit(1);
    }
    if (fastnc_options.threshold_name > 1.0) {
        print_help();
        fprintf(stderr,"\n%s: error: The threshold for deletion cannot be greater than 1.0\n", argv[0]);
        exit(1);
    }
	if (fastnc_options.number_name < 0) {
        print_help();
        fprintf(stderr,"\n%s: error: Number of iterations cannot be less than 0\n", argv[0]);
        exit(1);
    }
	
	// Check that the adj_table file exists
    std::ifstream adj_table;
    adj_table.open(fastnc_options.adj_table_filename);
    if (!adj_table.good()) {
        print_help();
        fprintf(stderr, "\n%s: error: adj_table %s does not exist\n", argv[0], fastnc_options.adj_table_filename.c_str());
        exit(1);
    }


    return fastnc_options;
}
