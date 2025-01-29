#ifndef VECTOR3D_H
#define VECTOR3D_H

class Vector3D
{
  public:
    float x, y, z;
    Vector3D(){x=y=z=0;}
    Vector3D(float a, float b, float c){x=a, y=b, z=c;}
    ~Vector3D() {}

    float GetMagnitude()
    {
      return sqrt(x*x + y*y + z*z);
    }

    Vector3D operator*(float num) const
    {
      return Vector3D(x*num, y*num, z*num);
    }

    Vector3D operator+(const Vector3D &vec) const
    {
      return Vector3D(x + vec.x, y + vec.y, z + vec.z);
    }

    Vector3D operator-(const Vector3D &vec) const
    {
      return Vector3D(x - vec.x, y - vec.y, z - vec.z);
    }

    Vector3D &operator=(const Vector3D &vec)
    {
      x = vec.x;
      y = vec.y;
      z = vec.z;
      return *this;
    }

    void NormalizeVector3D(void)
    {
      float mag = sqrt(x*x + y*y + z*z);
      x /= mag;
      y /= mag;
      z /= mag;
    }

    float DotVector3D(const Vector3D &vec) const
    {
      return x*vec.x + y*vec.y + z*vec.z;
    }

    Vector3D CrossVector3D(const Vector3D &vec) const
    {
      return Vector3D(y*vec.z - z*vec.y,
          z*vec.x - x*vec.z,
          x*vec.y - y*vec.x);
    }

    float AngleBetweenVector3Ds(Vector3D& vec)
    {
      return ( acos( DotVector3D(vec) / (GetMagnitude() * vec.GetMagnitude())) * (180./M_PI));
    }

    void Print(void)
    {
      std::cout << "v.x: " << x << ", v.y: " << y << ", v.z: " << z << std::endl;
    }
};

#endif
