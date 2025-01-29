#ifndef RAYINTERSECT_H
#define RAYINTERSECT_H

#include <iostream>

#include "Vector3D.h"

class Triangle
{
  public:
    Vector3D vertex0;
    Vector3D vertex1;
    Vector3D vertex2;
    Triangle() {}
    Triangle(Vector3D v0, Vector3D v1, Vector3D v2)
    {
      vertex0 = v0, vertex1 = v1, vertex2 = v2;
    }

    void Print(void)
    {
      std::cout << "Triangle information:" << std::endl;
      std::cout << "1st vertex: (" << vertex0.x << ", " << vertex0.y << ", " << vertex0.z << ")." << std::endl;
      std::cout << "2nd vertex: (" << vertex1.x << ", " << vertex1.y << ", " << vertex1.z << ")." << std::endl;
      std::cout << "3rd vertex: (" << vertex2.x << ", " << vertex2.y << ", " << vertex2.z << ")." << std::endl;
      std::cout << "End of triangle information." << std::endl;
    }
};

class Box
{
  public:
    Vector3D vertex[8];
    Box(Vector3D v0, Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4, Vector3D v5, Vector3D v6, Vector3D v7)
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


bool RayIntersectsTriangle(Vector3D rayOrigin,
    Vector3D rayVector,
    Triangle* inTriangle,
    Vector3D& outIntersectionPoint);

bool RayIntersectsCube(Vector3D rayOrigin,
    Vector3D rayVector,
    Box* inBox);

#endif

