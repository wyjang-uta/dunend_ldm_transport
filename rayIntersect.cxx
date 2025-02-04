#include "rayIntersect.h"
#include "vector3D.h"

bool RayIntersectsTriangle(vector3D rayOrigin,
    vector3D rayVector,
    Triangle* inTriangle,
    vector3D& outIntersectionPoint)
{
  const float EPSILON = 1e-7;
  vector3D vertex0 = inTriangle->vertex0;
  vector3D vertex1 = inTriangle->vertex1;
  vector3D vertex2 = inTriangle->vertex2;
  vector3D edge1, edge2, h, s, q;
  float a, f, u, v;
  edge1 = vertex1 - vertex0;
  edge2 = vertex2 - vertex0;
  h = rayVector.cross(edge2);
  a = edge1.dot(h);
  if( a > -EPSILON && a < EPSILON)
    return false;   // This ray is parallel to this triangle.
  f = 1.0/a;
  s = rayOrigin - vertex0;
  u = f * s.dot(h);
  if( u < 0.0 || u > 1.0 )
    return false;
  q = s.cross(edge1);
  v = f * rayVector.dot(q);
  if( v < 0.0 || u + v > 1.0 )
    return false;
  // At this stage we can compute t to find out where the intersection point is on the line.
  float t = f * edge2.dot(q);
  if( t > EPSILON ) // ray intersection
  {
    outIntersectionPoint = rayOrigin + rayVector * t;
    return true;
  }
  else    // This means that there is a line intersection but not a ray intersection.
    return false;
}

bool RayIntersectsCube(vector3D rayOrigin,
    vector3D rayVector,
    Box* inBox)
{
  bool front, front1, front2;
  bool back, back1, back2;
  bool top, top1, top2;
  bool bottom, bottom1, bottom2;
  bool left, left1, left2;
  bool right, right1, right2;

  Triangle frontTri1, frontTri2;
  Triangle backTri1, backTri2;
  Triangle topTri1, topTri2;
  Triangle bottomTri1, bottomTri2;
  Triangle leftTri1, leftTri2;
  Triangle rightTri1, rightTri2;

  vector3D intersection;
  // front test
  frontTri1 = Triangle(inBox->vertex[0], inBox->vertex[1], inBox->vertex[2]);
  frontTri2 = Triangle(inBox->vertex[1], inBox->vertex[2], inBox->vertex[3]);
  front1 = RayIntersectsTriangle(rayOrigin, rayVector, &frontTri1, intersection);
  front2 = RayIntersectsTriangle(rayOrigin, rayVector, &frontTri2, intersection);
  front = front1 || front2;

  // back test
  backTri1 = Triangle(inBox->vertex[4], inBox->vertex[5], inBox->vertex[6]);
  backTri2 = Triangle(inBox->vertex[5], inBox->vertex[6], inBox->vertex[7]);
  back1 = RayIntersectsTriangle(rayOrigin, rayVector, &backTri1, intersection);
  back2 = RayIntersectsTriangle(rayOrigin, rayVector, &backTri2, intersection);
  back = back1 || back2;

  // top test
  topTri1 = Triangle(inBox->vertex[0], inBox->vertex[4], inBox->vertex[1]);
  topTri2 = Triangle(inBox->vertex[4], inBox->vertex[1], inBox->vertex[5]);
  top1 = RayIntersectsTriangle(rayOrigin, rayVector, &topTri1, intersection);
  top2 = RayIntersectsTriangle(rayOrigin, rayVector, &topTri2, intersection);
  top = top1 || top2;

  // bottom test
  bottomTri1 = Triangle(inBox->vertex[2], inBox->vertex[3], inBox->vertex[6]);
  bottomTri2 = Triangle(inBox->vertex[3], inBox->vertex[6], inBox->vertex[7]);
  bottom1 = RayIntersectsTriangle(rayOrigin, rayVector, &bottomTri1, intersection);
  bottom2 = RayIntersectsTriangle(rayOrigin, rayVector, &bottomTri2, intersection);
  bottom = bottom1 || bottom2;

  // left test
  leftTri1 = Triangle(inBox->vertex[0], inBox->vertex[4], inBox->vertex[2]);
  leftTri2 = Triangle(inBox->vertex[4], inBox->vertex[2], inBox->vertex[6]);
  left1 = RayIntersectsTriangle(rayOrigin, rayVector, &leftTri1, intersection);
  left2 = RayIntersectsTriangle(rayOrigin, rayVector, &leftTri2, intersection);
  left = left1 || left2;

  // right test
  rightTri1 = Triangle(inBox->vertex[1], inBox->vertex[5], inBox->vertex[3]);
  rightTri2 = Triangle(inBox->vertex[5], inBox->vertex[3], inBox->vertex[7]);
  right1 = RayIntersectsTriangle(rayOrigin, rayVector, &rightTri1, intersection);
  right2 = RayIntersectsTriangle(rayOrigin, rayVector, &rightTri2, intersection);
  right = right1 || right2;

  int trueCounter = 0;

  if( front ) trueCounter++;
  if( back ) trueCounter++;
  if( top ) trueCounter++;
  if( bottom ) trueCounter++;
  if( left ) trueCounter++;
  if( right ) trueCounter++;

  if( trueCounter == 2 ) return true;
  else if( trueCounter > 2 )
  {
    std::cout << "something strange? hit counter is greater than 2; trueCounter = " << trueCounter << std::endl;
    return true;
  }
  else return false;
}
