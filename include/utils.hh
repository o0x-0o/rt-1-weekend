#pragma once

const double pi = 3.1415926535897;
const double infinity = std::numeric_limits<double>::infinity();

inline double random_double();
inline double random_double(double, double);

inline double clamp(double, double, double);

inline double degrees_to_radians(double);

Vec3 random_in_hemisphere(const Vec3&);
