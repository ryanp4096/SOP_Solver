#ifndef MAIN_H
#define MAIN_H

#include "reader.hpp"
#include <iostream>
#include <unordered_map>

int main(int argc, char **argv);

void single(char *path);

void compare(char *path1, char *path2);

void compare_node(TraceReader &r1, TraceReader &r2);

#endif
