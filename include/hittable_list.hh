#pragma once

#include <memory>

#include "vec3.hh"
#include "hittable.hh"
#include "hit_record.hh"

class HittableList : public Hittable {
public:
  std::vector<std::shared_ptr<Hittable>> objs;

public:
  HittableList() {}
  HittableList(std::shared_ptr<Hittable>);

  void clear();
  void add(std::shared_ptr<Hittable>);

  virtual bool hit(const Ray&,
                   double,
                   double,
                   HitRecord&) const override;
};
