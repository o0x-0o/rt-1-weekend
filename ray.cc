#include "ray.hh"
#include "vec3.hh"

Point3 Ray::origin()    const { return this->orig; }
Vec3   Ray::direction() const { return this->dir;  }

Point3 Ray::look_at(double t) const {
  return this->orig + t*this->dir;
}
