#include "vec3.hh"
#include "utils.hh"


/*
* Coordinate getters
* ---
* Here a macro is used to reduce code redundancy.
* These functions just gets [x, y, z] individual components from a 3-dimensional vector.
*/
#define CORD_GETTER(xyz, index) \
  double Vec3::xyz() const { return this->xyz[index]; }

CORD_GETTER(x, 0);
CORD_GETTER(y, 1);
CORD_GETTER(z, 2);

#undef CORD_GETTER


/*
* Negative vector function
* ---
* Negative-signs a vector's coordinates,
* returning a new vector with negative coordinates.
*/
Vec3 operator-() const {
  // Macro to reduce code redundancy.
  // t(n) = negative vector component with index n.
  #define t(n) \
    - this->xyz[n]

  return Vec3(
    t(0), t(1), t(2)
  );

  #undef t
}


/*
* Indexing operators.
* ---
* Here macros were not used because I didn't wanted to
* compromise much of the code's cleaness and readability.
*/
double  operator[](int index) const { return this->xyz[index]; }
double& operator[](int index)       { return this->xyz[index]; }


/*
 * Vec3xVec3 addiction and subtraction functions.
*/
Vec3& operator+=(const Vec3 &other) {
  // Macro for sum the nth individual component from this vector
  // with the nth individual component from the other vector.
  #define sum(n) \
    this->xyz[n] += other.xyz[n];

  sum(0);
  sum(1);
  sum(2);

  #undef sum

  return *this;
}

Vec3& operator-=(const Vec3 &other) {
  // Same principle as above macro but for subtraction.
  #define minus(n) \
    this->xyz[n] -= other.xyz[n];

  minus(0);
  minus(1);
  minus(2);

  #undef minus

  return *this;
}


/*
 * Vec3xScalar multiplication and division operators.
*/
Vec3& operator*=(const double t) {
  // Macro for multiply nth individual component from this 
  // vector by scalar t.
  #define o(n) \
    this->xyz[n] *= t;

  o(0);
  o(1);
  o(2);

  #undef o 

  return *this;
}

Vec3& operator/=(const double t) {
  #define o(n) \
    this->xyz[n] /= t;

  o(0);
  o(1);
  o(2);

  #undef o

  return *this;
}


/*
 * Util functions
 * ---
 * Here we got functions for 3-dimensional a vector's length
 * and length squared, as well as a numerical algorithm to
 * check for which the vector is near zero or not, and
 * random vector generators.
 */
double len() const { return std::sqrt(len_squared()); }
double len_squared() const {
  #define squared(n) \
    this->xyz[n] * this->xyz[n]

  return squared(0) + squared(1) + squared(2);

  #undef squared
}

bool near_zero() const {
  const auto s = 1e-8;

  #define o(n) \
    ( fabs(this->xyz[n]) < s )

  return o(0) && o(1) && o(2);

  #undef o
}

inline static Vec3 random() {
  return Vec3(
    random_double(),
    random_double(),
    random_double());
}

inline static Vec3 random(double min,
                          double max) {
  #define rmm \
    random_double(min, max)

  return Vec3(rmm, rmm, rmm);

  #undef rmm
}


/*
 * Vector Dot- and Cross-product funtions
 */
inline double vdot(const Vec3 &u,
                   const Vec3 &v) {
  #define o(n) \
    u.xyz[n] * v.xyz[n]

  return o(0) + o(1) + o(2);

  #undef o
}

inline Vec3 vcross(const Vec3 &u,
                   const Vec3 &v) {
  #define o(i, j) \
    u.xyz[i] * v.xyz[j] - u.xyz[j] * v.xyz[i]

  return Vec3(o(1, 2),
              o(2, 0),
              o(0, 1));
}

inline Vec3 vunit(Vec3 v) { return v / v.len(); }

inline Vec3 operator+(const Vec3 &u,
                      const Vec3 &v) {
  #define o(n) \
    u.xyz[n] + v.xyz[n]

  return Vec3(o(0), o(1), o(2));

  #undef o
}

inline Vec3 operator-(const Vec3 &u,
                      const Vec3 &v) {
  #define o(n) \
    u.xyz[n] - v.xyz[n]

  return Vec3(o(0), o(1), o(2));

  #undef o
}

inline Vec3 operator*(const Vec3 &u,
                      const Vec3 &v) {
  #define o(n) \
    u.xyz[n] * v.xyz[n]

  return Vec3(o(0), o(1), o(2));

  #undef o
}

inline Vec3 operator/(const Vec3 &u,
                      const Vec3 &v) {
  #define o(n) \
    u.xyz[n] / v.xyz[n] 

  return Vec3(o(0), o(1), o(2));

  #undef o
}


inline Vec3 operator*(double t,
                      const Vec3 &v) {
  #define o(n) \
    t*v.xyz[n]

  return Vec3(o(0), o(1), o(2));

  #undef o
}

inline Vec3 operator*(const Vec3 &v,
                      double t) {
  #define o(n) \
    v.xyz[n]*t 

  return Vec3(o(0), o(1), o(2));

  #undef o
}


inline Vec3 operator/(double t,
                      const Vec3 &v) {
  #define o(n) \
    (1/t) * v.xyz[n]

  return Vec3(o(0), o(1), o(2));

  #undef o
}

inline Vec3 operator/(const Vec3 &v,
                      double t) {
  #define o(n) \
    v.xyz[n] * (1/t)

  return Vec3(o(0), o(1), o(2));

  #undef o
}


inline std::ostream& operator<<(std::ostream &out,
                                const Vec3 &v) {
  #define o(n) \
    v.xyz[n] << ' '

  return out << o(0) << o(1) << o(2);
}


Vec3 refract(const Vec3 &uv,
             const Vec3 &n,
             double etai_over_etat) {
  double cos_theta = fmin(vdot(-uv, n),
                          1.0);

  Vec3 r_out_perp = etai_over_etat * (uv + cos_theta*n);
}
