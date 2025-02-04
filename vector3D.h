#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

class vector3D
{
  public:
    vector3D();
    vector3D(float a, float b, float c);
    vector3D(const vector3D& vec);
    ~vector3D();

    float getX() const { return x; }
    float getY() const { return y; }
    float getZ() const { return z; }

    void setX(float value) { x = value; }
    void setY(float value) { y = value; }
    void setZ(float value) { z = value; }

    float getMagnitude() const;
    vector3D  operator*(float num) const;
    vector3D  operator+(const vector3D  &vec) const;
    vector3D  operator-(const vector3D  &vec) const;
    vector3D& operator=(const vector3D  &vec);
    vector3D& operator=(const vector3D &&vec) noexcept;

    void normalize(void);
    float dot(const vector3D &vec) const;
    vector3D cross(const vector3D &vec) const;
    float angleBetween(vector3D& vec, bool inDegrees = false) const;
    void print(void);

  private:
    float x, y, z;

};

#endif
