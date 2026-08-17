#ifndef TREE_H
#define TREE_H

#include <vector>

struct TreeNode {
    int node;
    std::vector<TreeNode> children;

    unsigned long long enumerated_nodes;
    unsigned int recursive_children;
};

class TraceTree {
private:
    int max_depth;
    int node_limit;
    TreeNode root;
    TreeNode *current;
    int depth;
    std::vector<TreeNode *> stack;

    void print_tree_node(TreeNode &node, int indent);

public:
    TraceTree(int max_depth, int node_limit = 0)
        : max_depth(max_depth), root{.node = -1}, current(&root), depth(0), node_limit(node_limit) {}
    void initial_state(std::vector<int> &nodes);
    void add_node(int node);
    void edit_node(unsigned long long enumerated_nodes, unsigned int recursive_children);
    void finish_node();
    void print_subtree();
};

#endif