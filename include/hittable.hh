#pragma once

#include "vec3.hh"
#include "hit_record.hh"

class Hittable {
  public:
    virtual bool hit(const Ray&,
                     double t_min,
                     double t_max,
                     HitRecord&) const;
}
