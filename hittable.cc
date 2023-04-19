#include "hittable.hh"

#include "vec3.hh"
#include "hit_record.hh"

virtual bool Hittable::hit(const Ray &r,
                           double min,
                           double max,
                           HitRecord &rec) const = 0;
