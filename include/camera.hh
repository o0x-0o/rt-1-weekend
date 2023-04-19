#pragma once

#include "vec3.hh"

class Camera {
private:
  Point3 origin, lower_left_corner;
  Vec3 horizontal, vertical;

public:
  Camera(Point3,
         Point3,
         Vec3,
         double,
         double);

  Ray get_ray(double,
              double) const;
};
