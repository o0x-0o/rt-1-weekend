#pragma once

#include "ray.hh"
#include "vec3.hh"
#include "hit_record.hh"

class Material {
public:
  virtual bool scatter(const Ray&,
                       const HitRecord&,
                       Color&,
                       Ray&) const = 0;
};
