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
    TraceReader r1(path1);
    TraceReader r2(path2);
    r1.read_header();
    r2.read_header();

    while (!r1.file.eof() && !r2.file.eof()) {
        r1.parse_initial_state();
        r2.parse_initial_state();
        if (r1.current_path != r2.current_path) {
            std::cout << "!! different initial states" << std::endl;
            r1.parse_node();
            r2.parse_node();

        } else {
            compare_node(r1, r2);
        }

        r1.file.peek();
        r2.file.peek();
    }
    if (!r1.file.eof()) std::cout << "!! r1 has remaining workload" << std::endl;
    if (!r2.file.eof()) std::cout << "!! r2 has remaining workload" << std::endl;

    std::cout << "r1 enumerated nodes: " << r1.enumerated_nodes << std::endl;
    std::cout << "r1 ready nodes: " << r1.ready_nodes << std::endl;
    std::cout << "r1 recursive nodes: " << r1.recursive_nodes << std::endl;

    std::cout << "r2 enumerated nodes: " << r2.enumerated_nodes << std::endl;
    std::cout << "r2 ready nodes: " << r2.ready_nodes << std::endl;
    std::cout << "r2 recursive nodes: " << r2.recursive_nodes << std::endl;

    std::cout << "r1 only times: " << r1_only_times << std::endl;
    std::cout << "r1 only nodes: " << r1_only_nodes << std::endl;
    std::cout << "r2 only times: " << r2_only_times << std::endl;
    std::cout << "r2 only nodes: " << r2_only_nodes << std::endl;

    for (int i = 0; i < 11; i++) {
        if (r1_only_times_by_type[i] == 0 && r2_only_times_by_type[i] == 0) continue;
        std::cout << "PRUNED BY CODE " << i << std::endl;
        std::cout << "r1 only times: " << r1_only_times_by_type[i] << std::endl;
        std::cout << "r1 only nodes: " << r1_only_nodes_by_type[i] << std::endl;
        std::cout << "r2 only times: " << r2_only_times_by_type[i] << std::endl;
        std::cout << "r2 only nodes: " << r2_only_nodes_by_type[i] << std::endl;
    }

}

void compare_node(TraceReader &r1, TraceReader &r2) {
    r1.recursive_nodes++;
    r2.recursive_nodes++;
    bool x1, x2;

    std::unordered_map<int, TraceCode> differences;

    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCheck node1;
        NodeCheck node2;
        x1 = r1.parse_node_check(&node1);
        x2 = r2.parse_node_check(&node2);
        if (!x1 || !x2) break;
        r1.enumerated_nodes++;
        r2.enumerated_nodes++;
        if (node1.node != node2.node) {
            std::cout << "!! different items in check list" << std::endl;
            continue;
        }
        if (node1.action != node2.action) {
            // std::cout << "different actions for node " << node1.node << "" << std::endl;
            if (node1.action == TRACE_NO_PRUNE || node2.action == TRACE_NO_PRUNE) {
                if (node1.action == TRACE_NO_PRUNE) {
                    std::cout << "node " << node1.node << ": node1 not pruned | node2 pruned by " << node2.action << "" << std::endl;
                    differences.insert({node1.node, node2.action});
                } else {
                    std::cout << "node " << node1.node << ": node2 not pruned | node1 pruned by " << node1.action << "" << std::endl;
                    differences.insert({node2.node, node1.action});
                }
            }
        }
    }
    if (x1) std::cout << "!! r1 was not at end of list" << std::endl;
    if (x2) std::cout << "!! r2 was not at end of list" << std::endl;


    x1 = true;
    x2 = true;
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCall node1;
        NodeCall node2;
        while (x1) {
            x1 = r1.parse_node_call(&node1);
            if (x1) r1.ready_nodes++;
            if (x1 && differences.find(node1.node) != differences.end()) {
                unsigned long long en = 0;
                if (node1.action == TRACE_ENUMERATE) {
                    r1.current_path.push_back(node1.node);
                    en = r1.parse_node();
                    r1.current_path.pop_back();
                }
                std::cout << "parsing node " << node1.node << " only in r1 (" << en << " nodes)" << std::endl;
                r1_only_times++;
                r1_only_nodes += en;
                TraceCode action = differences[node1.node];
                r1_only_times_by_type[action]++;
                r1_only_nodes_by_type[action] += en;

            } else break;
        }
        while (x2) {
            x2 = r2.parse_node_call(&node2);
            if (x2) r2.ready_nodes++;
            if (x2 && differences.find(node2.node) != differences.end()) {
                unsigned long long en = 0;
                if (node2.action == TRACE_ENUMERATE) {
                    r2.current_path.push_back(node2.node);
                    en = r2.parse_node();
                    r2.current_path.pop_back();
                }
                std::cout << "parsing node " << node2.node << " only in r2 (" << en << " nodes)" << std::endl;
                r2_only_times++;
                r2_only_nodes += en;
                TraceCode action = differences[node2.node];
                r2_only_times_by_type[action]++;
                r2_only_nodes_by_type[action] += en;

            } else break;
        }

        if (!x1 && !x2) break;

        if (!x2) {
            std::cout << "!! r2 missing non-different node\n" << std::endl;
            unsigned long long en = 0;
            if (node1.action == TRACE_ENUMERATE) {
                r1.current_path.push_back(node1.node);
                en = r1.parse_node();
                r1.current_path.pop_back();
            }
            std::cout << "!! parsing node " << node1.node << " only in r1 (" << en << " nodes)" << std::endl;
            continue;
        } else if (!x1) {
            std::cout << "!! r1 missing non-different node\n" << std::endl;
            unsigned long long en = 0;
            if (node2.action == TRACE_ENUMERATE) {
                r2.current_path.push_back(node2.node);
                en = r2.parse_node();
                r2.current_path.pop_back();
            }
            std::cout << "!! parsing node " << node2.node << " only in r2 (" << en << " nodes)" << std::endl;
            continue;
        }

        if (node1.node != node2.node) {
            std::cout << "!! r1 and r2 have different nodes (different order?)" << std::endl;

            unsigned long long en = 0;
            if (node1.action == TRACE_ENUMERATE) {
                r1.current_path.push_back(node1.node);
                en = r1.parse_node();
                r1.current_path.pop_back();
            }
            std::cout << "!! parsing node " << node1.node << " only in r1 (" << en << " nodes)" << std::endl;

            en = 0;
            if (node2.action == TRACE_ENUMERATE) {
                r2.current_path.push_back(node2.node);
                en = r2.parse_node();
                r2.current_path.pop_back();
            }
            std::cout << "!! parsing node " << node2.node << " only in r2 (" << en << " nodes)" << std::endl;

            continue;
        }

        if (node1.action != TRACE_ENUMERATE || node2.action != TRACE_ENUMERATE) {
            if (node1.action != TRACE_ENUMERATE && node2.action != TRACE_ENUMERATE) {
                // both cancelled

            } else if (node1.action != TRACE_ENUMERATE) {
                // node1 cancelled
                unsigned long long en = 0;
                if (node2.action == TRACE_ENUMERATE) {
                    r2.current_path.push_back(node2.node);
                    en = r2.parse_node();
                    r2.current_path.pop_back();
                }
                std::cout << "node 1 cancelled by " << node1.action << ", parsing node " << node2.node << " only in r2 (" << en << " nodes)" << std::endl;
                r2_only_times++;
                r2_only_nodes += en;
                r2_only_times_by_type[node1.action]++;
                r2_only_nodes_by_type[node1.action] += en;

            } else if (node2.action != TRACE_ENUMERATE) {
                // node2 cancelled
                unsigned long long en = 0;
                if (node1.action == TRACE_ENUMERATE) {
                    r1.current_path.push_back(node1.node);
                    en = r1.parse_node();
                    r1.current_path.pop_back();
                }
                std::cout << "node 2 cancelled by " << node2.action << ", parsing node " << node1.node << " only in r1 (" << en << " nodes)" << std::endl;
                r1_only_times++;
                r1_only_nodes += en;
                r1_only_times_by_type[node1.action]++;
                r1_only_nodes_by_type[node1.action] += en;

            }
            continue;
        }

        r1.current_path.push_back(node1.node);
        r2.current_path.push_back(node2.node);

        compare_node(r1, r2);

        r1.current_path.pop_back();
        r2.current_path.pop_back();
    }
}