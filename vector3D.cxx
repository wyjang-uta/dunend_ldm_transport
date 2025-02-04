#include <iostream>

#include "vector3D.h"

vector3D::vector3D() : x(0), y(0), z(0) {}

vector3D::vector3D(float a, float b, float c) : x(a), y(b), z(c) {}

vector3D::vector3D(const vector3D& vec) : x(vec.x), y(vec.y), z(vec.z) {}

vector3D::~vector3D() {};

float vector3D::getMagnitude() const {
    return sqrt( x*x + y*y + z*z );
}

vector3D vector3D::operator*(float num) const {
    return vector3D( num*x, num*y, num*z );
}

vector3D vector3D::operator+(const vector3D& vec) const {
    return vector3D( x + vec.getX(), y + vec.getY(), z + vec.getZ() );
}

vector3D vector3D::operator-(const vector3D& vec) const {
    return vector3D( x - vec.getX(), y - vec.getY() , z - vec.getZ() );
}

vector3D& vector3D::operator=(const vector3D &vec) {
    if( this != &vec ) {
        x = vec.getX();
        y = vec.getY();
        z = vec.getZ();
    }
    return *this;
} 

vector3D& vector3D::operator=(const vector3D &&vec) noexcept {
    if( this != &vec ) {
        x = vec.x;
        y = vec.y;
        z = vec.z;
    }
    return *this;
} 

void vector3D::normalize(void) {
    float mag = sqrt(x*x + y*y + z*z);
    x /= mag;
    y /= mag;
    z /= mag;
}

float vector3D::dot(const vector3D &vec) const {
    return x * vec.getX() + y * vec.getY() + z * vec.getZ();
}

vector3D vector3D::cross(const vector3D &vec) const {
    return vector3D( y * vec.getZ() - z * vec.getY(),
                     z * vec.getX() - x * vec.getZ(),
                     x * vec.getY() - y * vec.getX() );
}

float vector3D::angleBetween(vector3D& vec, bool inDegrees) const {
    float angle = acos( dot(vec) / (getMagnitude() * vec.getMagnitude() ));
    return inDegrees ? angle * (180./M_PI) : angle;
}

void vector3D::print(void) {
    std::cout << "vector3D(" << this->getX() << ", " << this->getY() << ", " << this->getZ() << ")\n";
}
