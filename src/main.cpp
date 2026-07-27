#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.hpp"
#include <stdio.h>
#include <chrono>
#include <string>
#include <algorithm>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 4) {
        cout << "Usage: ./sop_solver <instance_path> <thread_count> <config_path>" << endl;
        exit(1);
    }
    string instance_path = argv[1];
    int thread_count = atoi(argv[2]);
    string config_path = argv[3];

    solver s;
    setpriority(PRIO_PROCESS, 0, -20);

    if (argc > 4) {
        string trace_path = argv[4];
        s.enable_trace(trace_path);
    }

    Config config = parse_config(config_path);
    s.assign_parameter(config);

    s.solve(instance_path, thread_count);

    return 0;
}