#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

// STB Image library headers
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Mathematical constant pi
#define M_PI 3.14159265358979323846

// Random number generator setup
static std::default_random_engine engine;
static std::uniform_real_distribution<double> uniform(0, 1);

// Function to calculate square of a number
static inline double sqr(double x) { return x * x; }

// Vector class representing a 3D point
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }

    // Overloaded subscript operator for accessing coordinates
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    // Overloaded addition assignment operator
    Vector& operator+=(const Vector& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    // Function to calculate square of the norm of the vector
    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    // Function to normalize the vector
    void normalize() {
        double n = sqrt(norm2());
        coord[0] /= n;
        coord[1] /= n;
        coord[2] /= n;
    }

    // Get the normalized vector without modifying the original vector
    Vector getNormalized() {
        Vector result = *this;
        result.normalize();
        return result;
    }

    // Negation operator to get the negated vector
    Vector operator-() const {
        return Vector(-coord[0], -coord[1], -coord[2]);
    }

    double coord[3];
};

// Overloaded arithmetic operators for vector operations
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const Vector& a, double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

// Overloaded division operator for vector and scalar
Vector operator/(const Vector& a, double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

// Dot product of two vectors
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Cross product of two vectors
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1],
                  a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    Vector random_direction_local(cos(2 * M_PI * r1) * sqrt(1 - r2), sin(2 * M_PI * r1) * sqrt(1 - r2), sqrt(r2));
    Vector random(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
    Vector tangent_1 = cross(N, random);
    tangent_1.normalize();
    Vector tangent_2 = cross(tangent_1, N);
  
    return random_direction_local[2] * N + random_direction_local[0] * tangent_1 + random_direction_local[1] * tangent_2;
}

// Ray class representing a ray in 3D space
class Ray {
public:
    Vector O; // Origin
    Vector u; // Direction

    Ray(const Vector& O, const Vector& u) : O(O), u(u) {}
};

// Sphere class representing a sphere in 3D space
class Sphere {
public:
    Vector C;      // Center
    double R;      // Radius
    Vector albedo; // Surface color
    bool is_mirror; // Whether the sphere behaves like a mirror
    bool is_transparent; // Whether the sphere is transparent

    Sphere(const Vector& c, double r, const Vector& a, bool is_mirror = false, bool is_transp = false)
        : C(c), R(r), albedo(a), is_mirror(is_mirror), is_transparent(is_transp) {}

    // Function to check for intersection between the sphere and a ray
    bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const {
        double a = 1;
        double b = 2 * dot(d.u, d.O - C);
        double c = (d.O - C).norm2() - R * R;

        double delta = b * b - 4 * a * c;
        if (delta < 0) return false;
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        if (t2 < 0) return false;
        if (t1 > 0)
            t = t1;
        else
            t = t2;

        P = d.O + t * d.u;
        N = (P - C).getNormalized();
        return true;
    }
};

// Scene class representing a collection of spheres in 3D space
class Scene {
public:
    std::vector<Sphere> objects; // Collection of spheres
    Sphere *light;       //pointer to a spheric light source
    double intensity_light;      // Intensity of the light source

    Scene() {}

    // Function to add a sphere to the scene
    void addSphere(const Sphere& sphere) {
        objects.push_back(sphere);
    }

    // Function to check for intersection between the scene and a ray
    bool intersect(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
        bool has_inter = false;
        min_t = 1E99;
        for (int i = 0; i < objects.size(); ++i) {
            Vector localP, localN;
            double t;
            bool local_has_inter = objects[i].intersect(d, localP, localN, t);
            if (local_has_inter) {
                has_inter = true;
                if (t < min_t) {
                    min_t = t;
                    P = localP;
                    N = localN;
                    sphere_id = i;
                }
            }
        }
        return has_inter;
    }
};

// Function to calculate color for a pixel using ray tracing
Vector getColor(Ray &r, const Scene &s, int nbrebonds) {
    if (nbrebonds == 0) return Vector(0, 0, 0);

    Vector P, N;
    int sphere_id;
    double t;
    bool has_inter = s.intersect(r, P, N, sphere_id, t);

    Vector intensite_pix(0, 0, 0);
    if (has_inter) {

        // Handling transparency
        if (s.objects[sphere_id].is_transparent) {
            // Refraction
            double n1 = 1;
            double n2 = 1.3;
            Vector normale_pour_transparence(N);
            if (dot(r.u, N) > 0) {
                n1 = 1.3;
                n2 = 1;
                normale_pour_transparence = -N;
            }

            double radical = 1 - sqr(n1 / n2) * (1 - sqr(dot(normale_pour_transparence, r.u)));
            if (radical > 0) {
                Vector direction_refraction = (n1 / n2) * (r.u - dot(r.u, normale_pour_transparence) * normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
                Ray rayon_refracte(P - 0.01 * normale_pour_transparence, direction_refraction);
                intensite_pix = getColor(rayon_refracte, s, nbrebonds - 1);
            }
        }
        // Handling mirrors
        else if (s.objects[sphere_id].is_mirror) {
            Vector direction_mirroir = r.u - 2 * dot(N, r.u) * N;
            Ray rayon_mirroir(P + 0.01 * N, direction_mirroir);
            intensite_pix = getColor(rayon_mirroir, s, nbrebonds - 1);
        }
        // Handling direct illumination
        else {
            // Ray ray_light(P + 0.01 * N, (s.position_light - P).getNormalized());
            // Vector P_light, N_light;
            // int sphere_id_light;
            // double t_light;
            // bool has_inter_light = s.intersect(ray_light, P_light, N_light, sphere_id_light, t_light);
            // double d_light2 = (s.position_light - P).norm2();
            // if (has_inter_light && t_light * t_light < d_light2) {
            //     intensite_pix = Vector(0, 0, 0);
            // }
            // else {
            //     intensite_pix = s.objects[sphere_id].albedo / M_PI * (s.intensity_light * std::max(0., dot((s.position_light - P).getNormalized(), N)) / (s.position_light - P).norm2());
            // }
        
        Vector axe_OP = (P-s.light->C).getNormalized();
        Vector random_direction = random_cos((P-s.light->C).getNormalized());
        Vector random_point = random_direction* s.light->R + s.light->C;
        Vector wi = (random_point - P).getNormalized();
        double d_light_squared = (random_point - P).norm2();
        Vector Np = random_direction;

        Ray ray_light(P + 0.01 * N, wi);
        Vector P_light, N_light;
        int sphere_id_light;
        double t_light;
        bool has_inter_light = s.intersect(ray_light, P_light, N_light, sphere_id_light, t_light);

    
        if (has_inter_light && t_light * t_light < 0.99 * d_light_squared) {
            intensite_pix = Vector(0, 0, 0);
        }
        else {
        
        intensite_pix = (s.intensity_light / (4*M_PI*d_light_squared) * std::max(0., dot(N, wi)) /dot(axe_OP,random_direction)) *s.objects[sphere_id].albedo;
        }
        // Handling indirect illumination
     
        Vector random_dir = random_cos(N);
        Ray random_rayon(P + 0.001 * N, random_dir);
        intensite_pix += getColor(random_rayon, s, nbrebonds - 1) * s.objects[sphere_id].albedo;
    }
}
    return intensite_pix;
}

int main() {
    int W = 512;
    int H = 512;
    const int number_of_rays = 8;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));

    Scene s;

    Sphere sphere_light(Vector(-5, 30, 40), 15, Vector(1.,1.,1.));

    Sphere s1(Vector(0, 0, -55), 10., Vector(1, 0, 0));
    Sphere s1b(Vector(-15, 0, -35), 10., Vector(1, 1, 0.2),false,true);
    Sphere s1c(Vector(15, 0, -75), 10., Vector(1, 0, 1),true);
    Sphere s2(Vector(0, -1000, 0), 960., Vector(0.0, 0.4, 0.14));
    Sphere s3(Vector(0, 1000, 0), 960., Vector(0.2, 0.2, 0.9));
    Sphere s4(Vector(-1000, 0, 0), 965., Vector(0.0, 0.2, 0.9));
    Sphere s5(Vector(1000, 0, 0), 965., Vector(0.9, 0.5, 0.7));
    Sphere s6(Vector(0, 0, -1000), 900., Vector(0.2, 0.1, 0.1));

    s.addSphere(sphere_light);

    s.addSphere(s1);
    s.addSphere(s1b);
    s.addSphere(s1c);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    s.light = &sphere_light;
    s.intensity_light = 5E9;
    Vector camera(0, 0 , 0);
    double focus_distance = 55;
    std::vector<unsigned char> image(W * H * 3, 0);

    int objectId;
    double best_t;

#pragma omp parallel for schedule (dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            Vector color(0., 0., 0.);
            for (int k = 0; k < number_of_rays; k++) {
               //methode Box Muller 
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = sqrt(-2*log(r1));
                double dx = R*cos(2*M_PI*r2);
                double dy = R*sin(2*M_PI*r2);

                double dx_aperture = (uniform(engine) - 0.5)*5.;
                double dy_aperture = (uniform(engine) - 0.5)*5.;

                Vector u(j - W / 2. + 0.5 + dx, -i + H / 2. - 0.5 + dy, -d);
                u.normalize();

                Vector destination = camera + focus_distance * u;
                Vector new_origin = camera + Vector(dx_aperture, dy_aperture, 0);
                Ray r(new_origin, (destination - new_origin).getNormalized());
                color += getColor(r, s, 5) / number_of_rays;
            } 
            // Apply gamma correction and store the color values in the image buffer
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2))); // BLUE
        }
    }

    // Write the rendered image to file
    stbi_write_png("image_4_3.png", W, H, 3, &image[0], 0);

    return 0;
}
