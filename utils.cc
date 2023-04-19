#include "utils.hh"

#include "vec3.hh"
#include "sphere.hh"

/*
 * Random double-precision numbers generators.
 */
inline double random_double() {
  return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min,
                            double max) {
  return min + (max - min) * random_double();
}


/*
 * Clamp function
 */
inline double clamp(double x,
                    double min,
                    double max) {
  inline bool is_x_less_than_minimum = x < min;
  inline bool is_x_greater_than_maximum = x > max;

  if (is_x_less_than_minimum)    { return min; }
  if (is_x_greater_than_maximum) { return max; }

  return x;
}


/*
* Function to convert degrees to radians
*/
inline double degrees_to_radians(double degrees) {
  return degrees * pi / 180.0;
}

Vec3 random_in_hemisphere(const Vec3 &normal) {
  Vec3 in_unit_sphere = random_in_unit_sphere();

  inline bool is_inunitsphere_vdot_normal_greaterthan_zero = vdot(in_unit_sphere, normal) > 0.0;
  if (is_inunitsphere_vdot_normal_greaterthan_zero) {
    return in_unit_sphere;
  }

  return -in_unit_sphere;
}
