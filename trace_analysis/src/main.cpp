#include "main.hpp"
#include <cstring>

static unsigned long long r1_only_times = 0;
static unsigned long long r1_only_nodes = 0;
static unsigned long long r2_only_times = 0;
static unsigned long long r2_only_nodes = 0;
static std::vector<unsigned long long> r1_only_times_by_type(11);
static std::vector<unsigned long long> r1_only_nodes_by_type(11);
static std::vector<unsigned long long> r2_only_times_by_type(11);
static std::vector<unsigned long long> r2_only_nodes_by_type(11);

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "Usage:\n" << "Analyze trace: ./sop_trace <path>\n" << "Print trace: ./sop_trace <path> -d <max_depth>\n" << "Compare traces: ./sop_trace <path1> <path2>" << std::endl;
        return 1;

    } else if (argc == 2) {
        single(argv[1]);

    } else if (argc >= 3) {
        if (argc >= 4 && strcmp(argv[2], "-d") == 0) {
            if (argc >= 6 && strcmp(argv[4], "-l") == 0) {
                single(argv[1], atoi(argv[3]), atoi(argv[5]));
            } else {
                single(argv[1], atoi(argv[3]));
            }

        } else {
            compare(argv[1], argv[2]);
        }
    }
    return 0;
}


void single(char *path, int max_depth, int node_limit) {
    TraceReader reader(path, max_depth, node_limit);

    reader.parse();

    reader.print_results();
}

void compare(char *path1, char *path2) {
    TraceCompare compare(path1, path2);
    
    compare.parse();

    compare.print_results();
}