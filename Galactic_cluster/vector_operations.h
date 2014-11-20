#pragma once
#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
#include <vector>
#include <cmath>
#include <vec3.h>


using std::vector;

vector<vec3> operator*(vector<vec3> a, double k);
//{
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        a[i] = a[i]*k;
//    }
//    return a;
//}

// Function for adding two std::vector<vec3>, with test for equal length
vector<vec3> operator+(vector<vec3> a, vector<vec3> b);
//{
//    if (a.size() != b.size())
//    {
//        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
//        return a;
//    }
//    vector<vec3> c = vector<vec3>(a.size());
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        c[i] = a[i] + b[i];
//    }
//    return c;
//}

// Function for subtracting two std::vector<vec3>, with test for equal length
vector<vec3> operator-(vector<vec3> a, vector<vec3> b);
//{
//    if (a.size() != b.size())
//    {
//        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
//        return a;
//    }
//    vector<vec3> c = vector<vec3>(a.size());
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        c[i] = a[i] - b[i];
//    }
//    return c;
//}

#endif // VECTOR_OPERATIONS_H
