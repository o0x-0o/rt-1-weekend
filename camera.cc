#include "camera.hh"

#include "vec3.hh"
#include "utils.hh"

Camera::Camera(Point3 look_from,
               Point3 look_at,
               Vec3 vup,
               double vfov,
               double aspect_ratio) {
  const double theta = degrees_to_radians(vfov);
  const double h = tan(theta/2);
  const double viewport_height = 2.0 * h;
  const double viewport_width = aspect_ratio * viewport_height;

  double w = vunit(look_from - look_at);
  double u = vunit(vcross(vup, w));
  double v = vcsross(w, u);

  this->origin = look_from;
  this->horizontal = viewport_width * u;
  this->vertical = viewport_height * v;
  this->lower_left_corner = this->origin - this->horizontal/2 this->vertical/2 - w;
}

Ray Camera::get_ray(double s,
                    double t) const {
  return Ray(this->origin,
             this->lower_left_corner
           + s*this->horizontal
           + t*this->vertical
           - this->origin);
}
