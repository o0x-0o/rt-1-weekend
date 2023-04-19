#pragma once

#include "material.hh"

#include "ray.hh"
#include "vec3.hh"
#include "hit_record.hh"

class Metal : public Material {
public:
  Color albedo;
  double fuzz;

public:
  Metal(const Color &a,
        double f) : albedo(a),
                    fuzz(f < 1 ? f : 1) {}

  virtual bool scatter(const Ray&,
                       const HitRecord&,
                       Color&,
                       Ray&) const override;
};
