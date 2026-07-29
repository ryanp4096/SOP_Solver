#include "main.hpp"

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
        std::cout << "Usage:\n" << "Analyze trace: ./sop_trace <path>\n" << "Compare traces: ./sop_trace <path1> <path2>" << std::endl;
        return 1;
    } else if (argc == 2) {
        single(argv[1]);
    } else if (argc >= 3) {
        compare(argv[1], argv[2]);
    }
    return 0;
}


void single(char *path) {
    TraceReader reader(path);

    reader.parse();

    std::cout << reader.enumerated_nodes << std::endl;
    std::cout << reader.ready_nodes << std::endl;
    std::cout << reader.recursive_nodes << std::endl;
}

void compare(char *path1, char *path2) {
    TraceCompare compare(path1, path2);
    
    compare.parse();

    compare.print_results();
}