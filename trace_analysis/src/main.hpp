#ifndef MAIN_H
#define MAIN_H

#include "reader.hpp"
#include "compare.hpp"
#include <iostream>

int main(int argc, char **argv);

void single(char *path, int max_depth = 0, int node_limit = 0);

void compare(char *path1, char *path2);

#endif
