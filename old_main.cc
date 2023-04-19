#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <limits>
#include <cstdlib>

#include "main.hh"

using std::sqrt;
using std::make_shared;
using std::shared_ptr;

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

inline double random_double() { return rand() / (RAND_MAX + 1.0); }
inline double random_double(double min, double max) {
	return min + (max-min)*random_double();
}

inline double clamp(double, double, double);
class Material;
class Hittable;
class Ray;

class Vec3 {
	public:
		Vec3() : e{0,0,0} {}
		Vec3(double e0,
		     double e1,
		     double e2) : e{e0,e1,e2} {}

		double x() const { return e[0]; }
		double y() const { return e[1]; }
		double z() const { return e[2]; }

		Vec3 operator-() const { return Vec3(-e[0], -e[1], e[2]); }
		double operator[](int i) const { return e[i]; }
		double& operator[](int i) { return e[i]; }

		Vec3& operator+=(const Vec3 &v) {
			e[0] += v.e[0];
			e[1] += v.e[1];
			e[2] += v.e[2];

			return *this;
		}

		Vec3& operator*=(const double t) {
			e[0] *= t;
			e[1] *= t;
			e[2] *= t;

			return *this;
		}

		Vec3& operator/=(const double t) {
			return *this *= 1/t;
		}

		double len() const {
			return sqrt(len_squared());
		}

		double len_squared() const {
			return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
		}

		bool near_zero() const {
			const auto s = 1e-8;
			return (fabs(e[0]) < s)
			    && (fabs(e[1]) < s)
			    && (fabs(e[2]) < s);
		}

		inline static Vec3 random() {
			return Vec3(random_double(),
				    random_double(),
				    random_double());
		}

		inline static Vec3 random(double min, double max) {
#define rmm \
			random_double(min, max)

			return Vec3(rmm, rmm, rmm);

#undef rmm
		}

	public:
		double e[3];
};

using Point3 = Vec3;
using Color = Vec3;

/*Vec3 random_in_unit_sphere() {
	while (true) {
		auto p = Vec3::random(-1, 1);
		if (p.len_squared() >= 1) { continue; }

		return p;
	}
}

Vec3 random_vunit() { return vunit(random_in_unit_sphere()); }*/

Color ray_color(const Ray& r, const Hittable& world, int depth);

inline std::ostream& operator<<(std::ostream &out, const Vec3 &v) {
	return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline Vec3 operator+(const Vec3& u, const Vec3& v) {
	return Vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline Vec3 operator-(const Vec3& u, const Vec3& v) {
	return Vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline Vec3 operator*(const Vec3& u, const Vec3& v) {
	return Vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline Vec3 operator*(double t, const Vec3& v) {
	return Vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline Vec3 operator*(const Vec3& v, double t) {
	return t * v;
}

inline Vec3 operator/(Vec3 v, double t) {
	return (1/t) * v;
}

inline double vdot(const Vec3& u, const Vec3& v) {
	return u.e[0] * v.e[0]
	     + u.e[1] * v.e[1]
	     + u.e[2] * v.e[2];
}

inline Vec3 vcross(const Vec3& u, const Vec3& v) {
	return Vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
		    u.e[2] * v.e[0] - u.e[0] * v.e[2],
		    u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline Vec3 vunit(Vec3 v) {
	return v / v.len();
}

Vec3 reflect(const Vec3& v, const Vec3& n) {
	return v - 2*vdot(v, n)*n;
}

Vec3 random_in_unit_sphere() {
        while (true) {
                auto p = Vec3::random(-1, 1);     
                if (p.len_squared() >= 1) { continue; }
                              
                return p;
        }                                    
}                                                      
                                                  
Vec3 random_vunit() { return vunit(random_in_unit_sphere()); }

void write_color(std::ostream& out, Color pixel, int samples_per_pixel) {
	auto r = pixel.x();
	auto g = pixel.y();
	auto b = pixel.z();

	auto scale = 1.0 / samples_per_pixel;
/*	r *= scale;
	g *= scale;
	b *= scale;*/
	r = sqrt(scale * r);
	g = sqrt(scale * g);
	b = sqrt(scale * b);

#define t(v) \
	static_cast<int>(256 * clamp(v, 0.0, 0.999))

	out << t(r) << ' ' << t(g) << ' ' << t(b) << '\n';

#undef t
}

class Ray {
	public:
		Ray() {}
		Ray(const Point3& origin, const Vec3& direction)
			: orig(origin), dir(direction) {}

		Point3 origin() const { return orig; }
		Vec3 direction() const { return dir; }

		Point3 at(double t) const {
			return orig + t*dir;
		}

	public:
		Point3 orig;
		Vec3 dir;
};

struct HitRecord {
	Point3 p;
	Vec3 normal;
	shared_ptr<Material> mat_ptr;
	double t;
	bool front_face;

	inline void set_face_normal(const Ray& r, const Vec3& outward_normal) {
		front_face = vdot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

class Hittable {
	public:
		virtual bool hit(const Ray& r, double t_min, double t_max,
				 HitRecord& rec) const = 0;
};

class Sphere : public Hittable {
	public:
		Sphere() {}
		Sphere(Point3 cen, double r,
		       shared_ptr<Material> m) : center(cen), radius(r), mat_ptr(m) {};

		virtual bool hit(const Ray& r, double t_min, double t_max,
				 HitRecord& rec) const override;
	
	public:
		Point3 center;
		double radius;
		shared_ptr<Material> mat_ptr;
};

bool Sphere::hit(const Ray& r, double t_min, double t_max,
		 HitRecord& rec) const {
	Vec3 oc = r.origin() - center;

	auto a = r.direction().len_squared();
	auto half_b = vdot(oc, r.direction());
	auto c = oc.len_squared() - radius*radius;
	
	auto discriminant = half_b*half_b - a*c;
	if (discriminant < 0) { return false; }

	auto sqrtd = sqrt(discriminant);

	auto root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root) {
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root) { return false; }
	}

	rec.t = root;
	rec.p = r.at(rec.t);
	Vec3 outward_normal = (rec.p - center) / radius;
	rec.set_face_normal(r, outward_normal);
	rec.mat_ptr = mat_ptr;

	return true;
}

class HittableList : public Hittable {
	public:
	HittableList() {}
	HittableList(shared_ptr<Hittable> obj) { add(obj); }
	
	void clear() { objs.clear(); }
	void add(shared_ptr<Hittable> obj) { objs.push_back(obj); }

	virtual bool hit(const Ray& r, double t_min, double t_max,
			 HitRecord& rec) const override;

	public:
	std::vector<shared_ptr<Hittable>> objs;
};

bool HittableList::hit(const Ray& r, double t_min, double t_max,
		       HitRecord& rec) const {
	HitRecord temp_rec;
	bool hit_anything = false;
	auto closest_so_far = t_max;

	for (const auto& obj: this->objs) {
		if (obj->hit(r, t_min, closest_so_far, temp_rec)) {
			hit_anything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}

	return hit_anything;
}

double hit_sphere(const Point3& center, double radius, const Ray& r) {
	Vec3 oc = r.origin() - center;
	auto a = r.direction().len_squared();
	auto half_b = vdot(oc, r.direction());
	auto c = oc.len_squared() - radius*radius;
	auto discriminant = half_b*half_b - a*c;

	if (discriminant < 0) {
		return -1.0;
	}

	return (-half_b - sqrt(discriminant) ) / a;
}

class Camera {
	public:
		Camera(Point3 look_from, Point3 look_at, Vec3 vup,
		       double vfov, double aspect_ratio) {
			auto theta = degrees_to_radians(vfov);
			auto h = tan(theta/2);
			auto viewport_height = 2.0 * h;
			auto viewport_width = aspect_ratio * viewport_height;
			
			auto w = vunit(look_from - look_at);
			auto u = vunit(vcross(vup, w));
			auto v = vcross(w, u);

			origin = look_from;
			horizontal = viewport_width * u;
			vertical = viewport_height * v;
			lower_left_corner = origin - horizontal/2 - vertical/2
					  - w;
		}

		Ray get_ray(double s, double t) const {
			return Ray(origin, lower_left_corner
					 + s*horizontal
					 + t*vertical
					 - origin);
		}

	private:
		Point3 origin;
		Point3 lower_left_corner;
		Vec3 horizontal;
		Vec3 vertical;
};

inline double clamp(double x, double min, double max) {
	if (x < min) { return min; }
	if (x > max) { return max; }

	return x;
}

Vec3 random_in_hemisphere(const Vec3& normal) {
	Vec3 in_unit_sphere = random_in_unit_sphere();

	if (vdot(in_unit_sphere, normal) > 0.0) { return  in_unit_sphere; }
	else 					{ return -in_unit_sphere; }
}

class Material {
	public:
		virtual bool scatter(const Ray& r_in, const HitRecord& rec,
				     Color& attenuation, Ray& scattered) const = 0;
};

class Metal : public Material {
	public:
		Metal(const Color& a, double f) : albedo(a),
						  fuzz(f < 1 ? f : 1) {}

		virtual bool scatter(const Ray& r_in, const HitRecord& rec,
				     Color& attenuation, Ray& scattered) const override {
			Vec3 reflected = reflect(vunit(r_in.direction()), rec.normal);
			scattered = Ray(rec.p, reflected + fuzz*random_in_unit_sphere());
			attenuation = albedo;

			return (vdot(scattered.direction(), rec.normal) > 0);
		}

	public:
		Color albedo;
		double fuzz;
};

Vec3 refract(const Vec3& uv, const Vec3& n, double etai_over_etat) {
	auto cos_theta = fmin(vdot(-uv, n), 1.0);
	Vec3 r_out_perp = etai_over_etat * (uv + cos_theta*n);
	Vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.len_squared())) * n;

	return r_out_perp + r_out_parallel;
}

class Dielectric : public Material {
	public:
		Dielectric(double index_of_refraction)
			: ir(index_of_refraction) {}

		virtual bool scatter(const Ray& r_in, const HitRecord& rec,
				     Color& attenuation, Ray& scattered) const override {
			attenuation = Color(1.0, 1.0, 1.0);
			double refrac_ratio = rec.front_face ? (1.0/ir) : ir;

			Vec3 unit_dir = vunit(r_in.direction());
			
			double cos_theta = fmin(vdot(-unit_dir, rec.normal), 1.0);
			double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

			bool cannot_refract = refrac_ratio * sin_theta > 1.0;
			Vec3 direction;

			if (cannot_refract
			 || reflectance(cos_theta, refrac_ratio) > random_double()) {
				direction = reflect(unit_dir, rec.normal);
			} else { refract(unit_dir, rec.normal, refrac_ratio); }

			scattered = Ray(rec.p, direction);

			return true;
		}

	public:
		double ir;
	
	private:
		static double reflectance(double cosine, double ref_idx) {
			auto r0 = (1-ref_idx) / (1+ref_idx);
			r0 = r0*r0;

			return r0 + (1-r0)*pow((1 - cosine), 5);
		}
};

Color ray_color(const Ray& r, const Hittable& world, int depth) {
	HitRecord rec;

	if (depth <= 0) { return Color(0, 0, 0); }

	if (world.hit(r, 0.001, infinity, rec)) {
		Ray scattered;
		Color attenuation;

		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
			return attenuation * ray_color(scattered, world, depth-1);
		}

		return Color(0, 0, 0);
	}

	Vec3 unit_dir = vunit(r.direction());
	auto t = 0.5 * (unit_dir.y() + 1.0);

	return (1.0-t) * Color(1.0, 1.0, 1.0) + t*Color(0.5, 0.7, 1.0);
}

class Lambertian : public Material {
	public:
		Lambertian(const Color& a) : albedo(a) {}

		virtual bool scatter(const Ray& r_in, const HitRecord& rec,
				     Color& attenuation, Ray& scattered) const override {
			auto scatter_direction = rec.normal + random_vunit();
			if (scatter_direction.near_zero()) {
				scatter_direction = rec.normal;
			}

			scattered = Ray(rec.p, scatter_direction);
			attenuation = albedo;

			return true;
		}
	
	public:
		Color albedo;
};

int main(void) {
	const auto aspect_ratio = 16.0 / 9.0;
	const int image_width = 1280;
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	const int samples_per_pixel = 100;
	const int max_depth = 50;

/*	HittableList world;

	auto material_ground = make_shared<Lambertian>(Color(0.8, 0.8, 0.0));
//	auto material_center = make_shared<Lambertian>(Color(0.7, 0.3, 0.3));
//	auto material_left   = make_shared<Metal>(Color(0.8, 0.8, 0.8), 0.3);
	auto material_center = make_shared<Dielectric>(1.5);
	auto material_left   = make_shared<Dielectric>(1.5);
	auto material_right  = make_shared<Metal>(Color(0.8, 0.6, 0.2), 1.0);

	world.add(make_shared<Sphere>(Point3(0, -100.5, -1.0), 100.0, material_ground));
	world.add(make_shared<Sphere>(Point3(0.0, 0.0, -1.0), 0.5, material_center));
	world.add(make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), 0.5, material_left));
	world.add(make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), -0.4, material_left));
	world.add(make_shared<Sphere>(Point3(1.0, 0.0, -1.0), 0.5, material_right));*/

//	auto R = cos(pi/4);
	HittableList world;

//	auto material_left  = make_shared<Lambertian>(Color(0, 0, 1));
//	auto material_right = make_shared<Lambertian>(Color(1, 0, 0));
	auto material_ground = make_shared<Lambertian>(Color(0.8, 0.8, 0.0));
	auto material_center = make_shared<Lambertian>(Color(0.1, 0.2, 0.5));
	auto material_left   = make_shared<Dielectric>(1.5);
	auto material_right  = make_shared<Metal>(Color(0.8, 0.6, 0.2), 0.0);

	world.add(make_shared<Sphere>(Point3(0.0, -100.5, -1.0), 100.0,
				      material_ground));
	world.add(make_shared<Sphere>(Point3(0.0, 0.0, -1.0), 0.5, material_center));
	world.add(make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), 0.5, material_left));
	world.add(make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), -0.45, material_left));
	world.add(make_shared<Sphere>(Point3(1.0, 0.0, -1.0), 0.5, material_right));

	Camera cam(Point3(-2,2,1), Point3(0,0,-1), Vec3(0,1,0), 90, aspect_ratio);

	auto viewport_height = 2.0;
	auto viewport_width = aspect_ratio * viewport_height;
	auto focal_length = 1.0;

	auto orig = Point3(0, 0, 0);
	auto horizontal = Vec3(viewport_width, 0, 0);
	auto vertical = Vec3(0, viewport_height, 0);
	auto lower_left_corner = orig - horizontal/2 - vertical/2 - Vec3(0, 0, focal_length);

	std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

	int i, j;
	for (j = image_height - 1; j >= 0; --j) {
		std::cerr << "\nScanlines remaining: " << j << " " << std::flush;
		for (i = 0; i < image_width; ++i) {
/*			auto r = double(i) / (image_width-1);
			auto g = double(j) / (image_height-1);
			auto b = 0.25;

			int ir = static_cast<int>(255.999 * r);
			int ig = static_cast<int>(255.999 * g);
			int ib = static_cast<int>(255.999 * b);

			std::cout << ir << ' ' << ig << ' ' << ib << '\n';

			auto u = double(i) / (image_width-1);
			auto v = double(j) / (image_height-1);

			Ray r(orig, lower_left_corner + u*horizontal + v*vertical - orig);
			Color pixel = ray_color(r, world);
			write_color(std::cout, pixel);*/

			Color pixel(0, 0, 0);
			for (int s = 0; s < samples_per_pixel; ++s) {
				auto u = (i + random_double()) / (image_width-1);
				auto v = (j + random_double()) / (image_height-1);
				Ray r = cam.get_ray(u, v);
				pixel += ray_color(r, world, max_depth);
			}

			write_color(std::cout, pixel, samples_per_pixel);
		}
	}

	std::cerr << "\nDone.\n";
	
	return 0;
}
