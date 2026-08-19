#include "compare.hpp"

TraceCompare::TraceCompare(const std::string& path1, const std::string& path2)
    : path1(path1), path2(path2), r1(path1), r2(path2) {}

void TraceCompare::print_results() {
    std::cout << std::endl;

    std::cout << "Results:" << std::endl;

    std::cout << "  r1 enumerated nodes: " << r1.enumerated_nodes << std::endl;
    std::cout << "  r1 ready nodes: " << r1.ready_nodes << std::endl;
    std::cout << "  r1 recursive nodes: " << r1.recursive_nodes << std::endl;

    std::cout << std::endl;

    std::cout << "  r2 enumerated nodes: " << r2.enumerated_nodes << std::endl;
    std::cout << "  r2 ready nodes: " << r2.ready_nodes << std::endl;
    std::cout << "  r2 recursive nodes: " << r2.recursive_nodes << std::endl;

    std::cout << std::endl;

    std::cout << "  shared nodes: " << shared_nodes << std::endl;
    std::cout << "  r1 only nodes: " << r1_only_nodes << " (times: " << r1_only_times << ")" << std::endl;
    std::cout << "  r2 only nodes: " << r2_only_nodes << " (times: " << r2_only_times << ")" << std::endl;
    if (r1_only_error_times > 0 || r2_only_error_times > 0) {
        std::cout << "  r1 only error nodes: " << r1_only_error_nodes << " (times: " << r1_only_error_times << ")" << std::endl;
        std::cout << "  r2 only error nodes: " << r2_only_error_nodes << " (times: " << r2_only_error_times << ")" << std::endl;
    }

    std::cout << std::endl;

    std::cout << "  r1 counted enumerated nodes: " << shared_nodes + r1_only_nodes + r1_only_error_nodes << " (expected: " << r1.enumerated_nodes << ")" << std::endl;
    std::cout << "  r2 counted enumerated nodes: " << shared_nodes + r2_only_nodes + r2_only_error_nodes << " (expected: " << r2.enumerated_nodes << ")" << std::endl;

    std::cout << std::endl;

    for (int i = 0; i < 11; i++) {
        if (r1_only_times_by_type[i] == 0 && r2_only_times_by_type[i] == 0) continue;
        std::cout << "  Pruned by " << TRACE_CODE_NAMES[i] << ":" << std::endl;
        std::cout << "    r1 only nodes: " << r1_only_nodes_by_type[i] << " (times: " << r1_only_times_by_type[i] << ")" << std::endl;
        std::cout << "    r2 only nodes: " << r2_only_nodes_by_type[i] << " (times: " << r2_only_times_by_type[i] << ")" << std::endl;
        std::cout << std::endl;
    }
}

void TraceCompare::print_tree() {
    if (node_results.size() == 0) return;
    if (!node_results.back().different && node_results.back().children.empty()) return;

    for (int i = 0; i < node_results.size() - 1; i++) {
        std::cout << node_results.at(i).node;
        std::cout << " ";
    }
    print_tree_node(node_results.back(), 0);
    std::cout << std::endl;
}

void TraceCompare::print_tree_node(ResultsTreeNode& node, int indent) {
    if (node.children.size() == 0) {
        std::cout << node.node << ":" << std::endl;
        for (int i = 0; i < indent + 1; i++) std::cout << "  ";
        std::cout << "- R1: " << node.r1 << std::endl;
        for (int i = 0; i < indent + 1; i++) std::cout << "  ";
        std::cout << "- R2: " << node.r2 << std::endl;

    } else if (node.children.size() == 1 && indent > 0) {
        std::cout << node.node << " ";
        print_tree_node(node.children.front(), indent);

    } else {
        std::cout << node.node << ":" << std::endl;
        for (ResultsTreeNode& n : node.children) {
            for (int i = 0; i < indent + 1; i++) std::cout << "  ";
            print_tree_node(n, indent + 1);
        }
    }
}

void TraceCompare::parse() {
    r1.read_header();
    r2.read_header();

    while (!r1.file.eof() && !r2.file.eof()) {
        r1.parse_initial_state();
        r2.parse_initial_state();
        if (r1.current_path != r2.current_path) {
            {
                node_results.clear();
                for (int i = 0; i < r1.current_path.size(); i++) node_results.push_back({.node = r1.current_path.at(i)});
                unsigned long long en = r1.parse_node().enumerated_nodes;

                r1_only_error_times++;
                r1_only_error_nodes += en;

                node_results.back().different = true;
                node_results.back().r1 = "Enumerated " + std::to_string(en) + " nodes";
                node_results.back().r2 = "!! Different initial state";
                print_tree();
            }
            {
                node_results.clear();
                for (int i = 0; i < r2.current_path.size(); i++) node_results.push_back({.node = r2.current_path.at(i)});
                unsigned long long en = r2.parse_node().enumerated_nodes;

                r2_only_error_times++;
                r2_only_error_nodes += en;

                node_results.back().different = true;
                node_results.back().r1 = "!! Different initial state";
                node_results.back().r2 = "Enumerated " + std::to_string(en) + " nodes";
                print_tree();
            }

        } else {
            node_results.clear();
            for (int i = 0; i < r1.current_path.size(); i++) node_results.push_back({.node = r1.current_path.at(i)});
            compare_node();
            print_tree();
        }

        r1.file.peek();
        r2.file.peek();
    }
    while (!r1.file.eof()) {
        node_results.clear();
        for (int i = 0; i < r1.current_path.size(); i++) node_results.push_back({.node = r1.current_path.at(i)});
        unsigned long long en = r1.parse_node().enumerated_nodes;

        r1_only_error_times++;
        r1_only_error_nodes += en;

        node_results.back().different = true;
        node_results.back().r1 = "Enumerated " + std::to_string(en) + " nodes";
        node_results.back().r2 = "!! End of file";
        print_tree();
    }
    while (!r2.file.eof()) {
        node_results.clear();
        for (int i = 0; i < r2.current_path.size(); i++) node_results.push_back({.node = r2.current_path.at(i)});
        unsigned long long en = r2.parse_node().enumerated_nodes;

        r2_only_error_times++;
        r2_only_error_nodes += en;

        node_results.back().different = true;
        node_results.back().r1 = "!! End of file";
        node_results.back().r2 = "Enumerated " + std::to_string(en) + " nodes";
        print_tree();
    }
}

void TraceCompare::compare_node() {
    r1.recursive_nodes++;
    r2.recursive_nodes++;

    std::vector<TraceCode> r1_actions = std::vector<TraceCode>(r1.instance_size);
    std::vector<TraceCode> r2_actions = std::vector<TraceCode>(r2.instance_size);

    int r1_en;
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCheck node;
        bool available = r1.parse_node_check(&node);
        if (!available) break;
        r1.enumerated_nodes++;
        r1_en++;
        r1_actions[node.node] = node.action;
    }
    int r2_en;
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCheck node;
        bool available = r2.parse_node_check(&node);
        if (!available) break;
        r2.enumerated_nodes++;
        r2_en++;
        r2_actions[node.node] = node.action;
    }
    if (r1_en != r2_en) {
        std::cout << "!! different number of enumerated nodes" << std::endl;
    }
    shared_nodes += r1_en;

    bool x1 = true;
    bool x2 = true;
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCall node1;
        NodeCall node2;

        /* Go through nodes only in r1 until a node in both r1 and r2 is reached (or end of list) */
        while (x1) {
            x1 = r1.parse_node_call(&node1);
            if (!x1) break;
            r1.ready_nodes++;
            if (r2_actions[node1.node] != TRACE_NO_PRUNE) {
                if (node1.action == TRACE_ENUMERATE) {
                    r1.current_path.push_back(node1.node);
                    unsigned long long en = r1.parse_node().enumerated_nodes;
                    r1.current_path.pop_back();

                    r1_only_times++;
                    r1_only_nodes += en;
                    TraceCode action = r2_actions[node1.node];
                    r1_only_times_by_type[action]++;
                    r1_only_nodes_by_type[action] += en;

                    node_results.back().children.push_back({
                        .node = node1.node,
                        .different = true,
                        .r1 = "Enumerated " + std::to_string(en) + " nodes",
                        .r2 = "Pruned by " + TRACE_CODE_NAMES[action]
                    });
                }
            } else break;
        }

        /* Go through nodes only in r2 until a node in both r1 and r2 is reached (or end of list) */
        while (x2) {
            x2 = r2.parse_node_call(&node2);
            if (!x2) break;
            r2.ready_nodes++;
            if (r1_actions[node2.node] != TRACE_NO_PRUNE) {
                if (node2.action == TRACE_ENUMERATE) {
                    r2.current_path.push_back(node2.node);
                    unsigned long long en = r2.parse_node().enumerated_nodes;
                    r2.current_path.pop_back();

                    r2_only_times++;
                    r2_only_nodes += en;
                    TraceCode action = r1_actions[node2.node];
                    r2_only_times_by_type[action]++;
                    r2_only_nodes_by_type[action] += en;

                    node_results.back().children.push_back({
                        .node = node2.node,
                        .different = true,
                        .r1 = "Pruned by " + TRACE_CODE_NAMES[action],
                        .r2 = "Enumerated " + std::to_string(en) + " nodes",
                    });
                }
            } else break;
        }

        /* both r1 and r2 are out of items, done */
        if (!x1 && !x2) break;

        if (!x2) {
            /* r2 has ended list. r1 still has items */
            if (node1.action == TRACE_ENUMERATE) {
                r1.current_path.push_back(node1.node);
                unsigned long long en = r1.parse_node().enumerated_nodes;
                r1.current_path.pop_back();

                r1_only_error_times++;
                r1_only_error_nodes += en;

                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Enumerated " + std::to_string(en) + " nodes",
                    .r2 = "!! Missing"
                });
            } else {
                r1_only_error_times++;
                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Cancelled by " + TRACE_CODE_NAMES[node1.action],
                    .r2 = "!! Missing"
                });
            }
            continue;

        } else if (!x1) {
            /* r1 has ended list. r2 still has items */
            if (node2.action == TRACE_ENUMERATE) {
                r2.current_path.push_back(node2.node);
                unsigned long long en = r2.parse_node().enumerated_nodes;
                r2.current_path.pop_back();

                r2_only_error_times++;
                r2_only_error_nodes += en;

                node_results.back().children.push_back({
                    .node = node2.node,
                    .different = true,
                    .r1 = "!! Missing",
                    .r2 = "Enumerated " + std::to_string(en) + " nodes"
                });
            } else {
                r2_only_error_times++;
                node_results.back().children.push_back({
                    .node = node2.node,
                    .different = true,
                    .r1 = "!! Missing",
                    .r2 = "Cancelled by " + TRACE_CODE_NAMES[node2.action]
                });
            }
            continue;
        }

        if (node1.node != node2.node) {
            /* r1 and r2 both have this item, but are misaligned (usually different order) */

            if (node1.action == TRACE_ENUMERATE) {
                r1.current_path.push_back(node1.node);
                unsigned long long en = r1.parse_node().enumerated_nodes;
                r1.current_path.pop_back();

                r1_only_error_times++;
                r1_only_error_nodes += en;

                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Enumerated " + std::to_string(en) + " nodes",
                    .r2 = "!! Misaligned"
                });
            } else {
                r1_only_error_times++;
                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Cancelled by " + TRACE_CODE_NAMES[node1.action],
                    .r2 = "!! Misaligned"
                });
            }

            if (node2.action == TRACE_ENUMERATE) {
                r1.current_path.push_back(node1.node);
                unsigned long long en = r1.parse_node().enumerated_nodes;
                r1.current_path.pop_back();

                r1_only_error_times++;
                r1_only_error_nodes += en;

                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Enumerated " + std::to_string(en) + " nodes",
                    .r2 = "!! Misaligned"
                });
            } else {
                r1_only_error_times++;
                node_results.back().children.push_back({
                    .node = node2.node,
                    .different = true,
                    .r1 = "!! Missing",
                    .r2 = "Cancelled by " + TRACE_CODE_NAMES[node2.action]
                });
            }

            continue;
        }

        /* Both r1 and r2 have the same node, but either one is cancelled */
        if (node1.action != TRACE_ENUMERATE || node2.action != TRACE_ENUMERATE) {
            if (node1.action != TRACE_ENUMERATE && node2.action != TRACE_ENUMERATE) {
                /* Both r1 and r2 cancelled the node */

            } else if (node2.action != TRACE_ENUMERATE) {
                /* r2 cancelled the node but r1 didn't */
                r1.current_path.push_back(node1.node);
                unsigned long long en = r1.parse_node().enumerated_nodes;
                r1.current_path.pop_back();

                r1_only_times++;
                r1_only_nodes += en;
                r1_only_times_by_type[node2.action]++;
                r1_only_nodes_by_type[node2.action] += en;

                node_results.back().children.push_back({
                    .node = node1.node,
                    .different = true,
                    .r1 = "Enumerated " + std::to_string(en) + " nodes",
                    .r2 = "Cancelled by " + TRACE_CODE_NAMES[node2.action]
                });

            } else if (node1.action != TRACE_ENUMERATE) {
                /* r1 cancelled the node but r2 didn't */
                r2.current_path.push_back(node2.node);
                unsigned long long en = r2.parse_node().enumerated_nodes;
                r2.current_path.pop_back();

                r2_only_times++;
                r2_only_nodes += en;
                r2_only_times_by_type[node1.action]++;
                r2_only_nodes_by_type[node1.action] += en;

                node_results.back().children.push_back({
                    .node = node2.node,
                    .different = true,
                    .r1 = "Cancelled by " + TRACE_CODE_NAMES[node1.action],
                    .r2 = "Enumerated " + std::to_string(en) + " nodes"
                });
            }

            continue;
        }

        /* Both r1 and r2 have the same node and both are enumerated. Recursively call compare_node to compare this shared node*/

        node_results.push_back({.node = node1.node});

        r1.current_path.push_back(node1.node);
        r2.current_path.push_back(node2.node);

        compare_node();

        r1.current_path.pop_back();
        r2.current_path.pop_back();

        /* If this node has differences, add it as a child of its parent path in results tree*/
        if (!node_results.back().children.empty()) {
            node_results.at(node_results.size() - 2).children.push_back(node_results.back());
        }
        node_results.pop_back();
    }
}