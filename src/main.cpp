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
    setpriority(PRIO_PROCESS, 0, -20);

    string instance_path = argv[1];
    int thread_count = atoi(argv[2]);
    string config_path = argv[3];

    solver s(instance_path, thread_count, config_path);

    if (argc > 4) {
        string trace_path = argv[4];
        s.enable_trace(trace_path);
    }

    s.solve();

    return 0;
}