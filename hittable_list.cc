#include "hittable_list.hh"

bool HittableList::hit(const Ray &r,
                       double min,
                       double max,
                       HitRecord &rec) const {
  HitRecord temp_rec;
  bool hit_anything = false;
  double closest_so_far = max;

  for (const auto &obj: this->objs) {
    inline bool has_hitted_something = obj->hit(r, min, closest_so_far, temp_rec);

    if (has_hitted_something) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }

  return hit_anything;
}
