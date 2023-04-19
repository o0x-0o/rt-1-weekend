#include "metal.hh"

#include "vec3.hh"
#include "ray.hh"
#include "hit_record.hh"

bool scatter(const Ray &r_in,
             const HitRecord &rec,
             Color &attenuation,
             Ray &scattered) const override {
  Vec3 reflected = reflect(vunit(r_in.direction()),
                           rec.normal);

  scattered = Ray(rec.p,
                  reflected + this->fuzz*random_in_unit_sphere());
  attenuation = this->albedo;

  return vdot(scattered.direction(),
              rec.normal) > 0;
}
