#pragma once

class Vec3 {
public:
  union {
    // Vector coordinates / 
    double xyz[3];
    double rgb[3]; };

public:
  Vec3() : xyz{0, 0, 0} {}
  Vec3(double x,
       double y,
       double z) : xyz{x, y, z} {}

  // Coordinates getters
  double x() const;
  double y() const;
  double z() const;

  // Negative Vector
  Vec3 operator-() const;

  // Indexing
  double  operator[](int) const;
  double& operator[](int);

  // Vec3xVec3 addiction and subtraction
  Vec3& operator+=(const Vec3&);
  Vec3& operator-=(const Vec3&);

  // Vec3xScalar multiplication and division
  Vec3& operator*=(const double);
  Vec3& operator/=(const double);

  // Util functions
  double len() const;
  double len_squared() const;

  bool near_zero() const;

  inline static Vec3 random();
  inline static Vec3 random(double min,
                            double max);
};


// Vector cross product and Vector dot product
inline double vdot(const Vec3&,
                   const Vec3&);
inline Vec3 vcross(const Vec3&,
                   const Vec3&);

// Vector unit
inline Vec3 vunit(Vec3);

// Vec3xVec3 operations
inline Vec3 operator+(const Vec3&,
                      const Vec3&);
inline Vec3 operator-(const Vec3&,
                      const Vec3&);
inline Vec3 operator*(const Vec3&,
                      const Vec3&);
inline Vec3 operator/(const Vec3&,
                      const Vec3&);

// Vec3xScalar multiplication and division
inline Vec3 operator*(double,
                      const Vec3&);
inline Vec3 operator*(const Vec3&,
                      double);

inline Vec3 operator/(double,
                      const Vec3&);
inline Vec3 operator/(const Vec3&,
                      double);


// Vec3 to String 
inline std::ostream& operator<<(std::ostream&,
                                const Vec3&);


Vec3 refract(const Vec3&,
             const Vec3&,
             double);


// Aliases
using Point3 = Vec3;
using Color  = Vec3;
