#include "tree.hpp"

#include <iostream>

void TraceTree::initial_state(std::vector<int> &nodes) {
    if (max_depth < 1) return;

    current = &root;
    depth = 0;
    stack.clear();

    for (int node : nodes) {
        TreeNode *n = nullptr;
        for (TreeNode &child : current->children) {
            if (child.node == node) n = &child;
        }
        if (n == nullptr) {
            current->children.push_back(TreeNode{.node = node});
            n = &current->children.back();
        }
        current = n;
        depth++;
        stack.push_back(current);
    }
}

void TraceTree::add_node(int node) {
    if (max_depth < 1) return;

    if (depth >= max_depth || current == nullptr) {
        current = nullptr;
        depth++;
        stack.push_back(nullptr);
        return;
    }

    current->children.push_back(TreeNode{.node = node});
    current = &current->children.back();
    depth++;
    stack.push_back(current);
}

void TraceTree::finish_node() {
    if (max_depth < 1) return;
    if (node_limit > 0 && current != nullptr && current->enumerated_nodes < node_limit) {
        current->children.clear();
    }

    stack.pop_back();
    current = stack.back();
    depth--;
}

void TraceTree::edit_node(unsigned long long enumerated_nodes, unsigned int recursive_children) {
    if (max_depth < 1) return;
    if (current == nullptr) return;
    current->enumerated_nodes = enumerated_nodes;
    current->recursive_children = recursive_children;
}

void TraceTree::print_subtree() {
    if (max_depth < 1) return;

    for (int i = 0; i < stack.size() - 1; i++) {
        std::cout << stack.at(i)->node;
        std::cout << " ";
    }
    print_tree_node(*current, 0);
    std::cout << std::endl;
}

void TraceTree::print_tree_node(TreeNode &node, int indent) {
    if (node.children.size() == 0) {
        std::cout << node.node << ": " << node.recursive_children << " children, " << node.enumerated_nodes << " nodes" << std::endl;

    } else if (node.children.size() == 1 && indent > 0) {
        std::cout << node.node << " ";
        print_tree_node(node.children.front(), indent);

    } else {
        std::cout << node.node << ":" << std::endl;
        for (TreeNode &n : node.children) {
            for (int i = 0; i < indent + 1; i++) std::cout << "  ";
            print_tree_node(n, indent + 1);
        }
    }
}