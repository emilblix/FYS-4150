#pragma once
#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
#include <vector>
#include <cmath>
#include <vec3.h>

using std::vector;

vector<vec3> operator*(vector<vec3> a, double k);
vector<vec3> operator+(vector<vec3> a, vector<vec3> b);
vector<vec3> operator-(vector<vec3> a, vector<vec3> b);

#endif // VECTOR_OPERATIONS_H
