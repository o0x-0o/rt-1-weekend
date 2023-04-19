#pragma once

#include "vec3.hh"

class Ray
{
  public:
    Point3 orig; // Ray origin in space
    Vec3    dir; // Ray direction in space

  public:
    Ray() {}
    Ray(const Point3& origin,
        const   Vec3& direction) : orig(origin),
                                   dir(direction) {}

    Point3 origin() const;
    Vec3 direction() const;

    Point3 look_at(double) const;
};
