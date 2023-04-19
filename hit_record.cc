#include "hit_record.hh"

#include "vec3.hh"

inline void set_face_normal(const Ray &r,
                            const Vec3 &outward_normal) {
  this.front_face = vdot(r.direction(),
                         outward_normal) < 0;

  normal = this.front_face ? outward_normal : -outward_normal;
}
