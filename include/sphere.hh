#pragma once

#include <memory>

#include "vec3.hh"
#include "hittable.hh"

class Sphere : public Hittable {
  public:
    Point3 center;
    double radius;
    std::shared_ptr<Material> material;

  public:
    Sphere() {}
    Sphere(Point3 c,
           double r,
           std::shared_ptr<Material> m) : center(c),
                                          radius(r),
                                          material(m) {};

    virtual bool hit(const Ray&,
                     double,
                     double,
                     HitRecord&) const override;
};


Vec3 random_in_unit_sphere();

double hit_sphere(const Point3&,
                  double,
                  const Ray&);
