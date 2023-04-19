#pragma once

#include <memory>

#include "vec3.hh"

struct HitRecord
{
  Point3 p;
  Vec3 normal;
  std::shared_ptr<Material> material;
  double t;
  bool front_face;

  inline void set_face_normal(const Ray&,
                              const Vec3&);
};
