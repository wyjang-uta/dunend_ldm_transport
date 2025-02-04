#ifndef RAYINTERSECT_H
#define RAYINTERSECT_H

#include <iostream>

#include "vector3D.h"

class Triangle
{
  public:
    vector3D vertex0;
    vector3D vertex1;
    vector3D vertex2;
    Triangle() {}
    Triangle(vector3D v0, vector3D v1, vector3D v2)
    {
      vertex0 = v0, vertex1 = v1, vertex2 = v2;
    }

    void Print(void)
    {
      std::cout << "Triangle information:" << std::endl;
      std::cout << "1st vertex: (" << vertex0.getX() << ", " << vertex0.getY() << ", " << vertex0.getZ() << ")." << std::endl;
      std::cout << "2nd vertex: (" << vertex1.getX() << ", " << vertex1.getY() << ", " << vertex1.getZ() << ")." << std::endl;
      std::cout << "3rd vertex: (" << vertex2.getX() << ", " << vertex2.getY() << ", " << vertex2.getZ() << ")." << std::endl;
      std::cout << "End of triangle information." << std::endl;
    }
};

class Box
{
  public:
    vector3D vertex[8];
    Box(vector3D v0, vector3D v1, vector3D v2, vector3D v3, vector3D v4, vector3D v5, vector3D v6, vector3D v7)
    {
      vertex[0] = v0;
      vertex[1] = v1;
      vertex[2] = v2;
      vertex[3] = v3;
      vertex[4] = v4;
      vertex[5] = v5;
      vertex[6] = v6;
      vertex[7] = v7;
    }
};


bool RayIntersectsTriangle(vector3D rayOrigin,
    vector3D rayVector,
    Triangle* inTriangle,
    vector3D& outIntersectionPoint);

bool RayIntersectsCube(vector3D rayOrigin,
    vector3D rayVector,
    Box* inBox);

#endif

