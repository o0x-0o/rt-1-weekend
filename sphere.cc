#include "sphere.hh"

#include "vec3.hh"
#include "hit_record.hh"

bool Sphere::hit(const Ray &r,
                 double min,
                 double max,
                 HitRecord &rec) {
  Vec3 oc = r.origin() - center;

  double a = r.direction().len_squared();
  double half_b = vdot(oc, r.direction());
  double c = oc.len_squared() - this->radius*this->radius;

  double delta = half_b*half_b - a*c;

  inline bool is_delta_less_than_zero = delta < 0;
  if (is_delta_less_than_zero) { return false; }

  double sqrtd = sqrt(delta);

  double root = (-half_b - sqrtd) / a;
  
  inline bool is_root_less_than_minimum = root < min;
  inline bool is_root_greater_than_maximum = root > max;
  if (is_root_less_than_minimum || is_root_greater_than_maximum) {
    root = (-half_b + sqrtd) / a;
    return false;
  }

  rec.t = root;
  rec.p = r.look_at(rec.t);
  Vec3 outward_normal = (rec.p - center) / radius;
  rec.set_face_normal(r, outward_normal);
  rec.material = this->material;

  return true;
}


Vec3 random_in_unit_sphere() {
  while (true) {
    auto p = Vec3::random(-1, 1);
    inline bool is_p_squared_greater_than_one = p.len_squared() >= 1;

    if (is_p_squared_greater_than_one) { continue; }

    return p;
  }
}


double hit_sphere(const Point3 &center,
                  double radius,
                  const Ray &r) {
  Vec3 oc = r.origin() - center;

  double a = r.direction().len_squared();
  double half_b = vdot(oc, r.direction());
  double c = oc.len_squared() - radius*radius;
  double delta = half_b*half_b - a*c;

  inline bool is_delta_less_than_zero = delta < 0;
  if (is_delta_less_than_zero) {
    return -1.0;
  }

  return (-half_b - sqrt(delta)) / a;
}
