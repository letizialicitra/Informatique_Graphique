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


// Object class representing a generic object in the scene
class Object{
public:
    Object(){};
    virtual bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const = 0;

    Vector albedo; // Surface color
    bool is_mirror; // Whether the object behaves like a mirror
    bool is_transparent; // Whether the object is transparent

};


// Sphere class representing a sphere in 3D space
class Sphere : public Object {
public:
    Vector C;      // Center
    double R;      // Radius
    
    Sphere(const Vector& c, double r, const Vector& a, bool mirror = false, bool transp = false)
        : C(c), R(r) {
        albedo = a;
        is_mirror = mirror;
        is_transparent = transp;
        };

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

// Triangle class representing a triangle in 3D space
class Triangle : public Object{
public:
    // Constructor for Triangle class
    Triangle(const Vector& A, const Vector& B, const Vector& C,const Vector& a, bool mirror = false, bool transp = false) : A(A), B(B), C(C) {
        albedo = a; // Assign surface color
        is_mirror = mirror; // Set whether it behaves like a mirror
        is_transparent = transp; // Set whether it is transparent
    };

    // Function to check for intersection between the triangle and a ray
    bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const {

        // Calculate the normal to the triangle
        N = cross(B-A, C-A).getNormalized();
        t = dot(C-d.O, N) / dot(d.u,N); // Calculate intersection parameter t
        if (t<0) return false; // If intersection is behind the ray, return false

        P = d.O + t*d.u; // Calculate the intersection point
        Vector u = B-A; // Vector from A to B
        Vector v = C-A; // Vector from A to C
        Vector w = P-A; // Vector from A to the intersection point
        double m11 = (B-A).norm2(); // Squared length of u
        double m12 = dot(u,v); // Dot product of u and v
        double m22 = v.norm2(); // Squared length of v
        double det_m = m11*m22 - m12*m12; // Determinant of the matrix formed by u and v

        double b11 = dot(w,u); // Dot product of w and u
        double b21 = dot(w,v); // Dot product of w and v
        double det_b = b11*m22 - b21*m12; // Determinant of the matrix formed by w and v
        double beta = det_b /det_m; // Coordinate barycenter w.r.t to B

        double g12 = b11; // Dot product of w and u
        double g22 = b21; // Dot product of w and v
        double det_g = m11*g22 - m12*g12; // Determinant of the matrix formed by w and u
        double gamma = det_g / det_m; // Coordinate barycenter w.r.t to C

        double alpha = 1 - beta - gamma; // Calculate the alpha coordinate
        if ( alpha < 0 || alpha > 1 ) return false; // Check if alpha is within range
        if (beta < 0 || beta > 1 ) return false; // Check if beta is within range
        if (gamma <0 || gamma>1 ) return false; // Check if gamma is within range

        return true; // Return true if the intersection point is within the triangle
    }
    Vector A, B ,C; // Vertices of the triangle
};


// Scene class representing a collection of spheres in 3D space
class Scene {
public:
    std::vector<const Object *> objects; // Collection of objects
    Sphere *light;       //pointer to a spherical light source
    double intensity_light;      // Intensity of the light source

    Scene() {}

    // Function to add a sphere to the scene
    void addSphere(const Sphere& sphere) {
        objects.push_back(&sphere);
    }

    // Function to add a triangle to the scene
    void addTriangle(const Triangle& triangle) {
        objects.push_back(&triangle);
    }

    // Function to check for intersection between the scene and a ray
    bool intersect(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
        bool has_inter = false;
        min_t = 1E99;
        for (int i = 0; i < objects.size(); ++i) {
            Vector localP, localN;
            double t;
            bool local_has_inter = objects[i]->intersect(d, localP, localN, t);
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
        if (s.objects[sphere_id]->is_transparent) {
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
        else if (s.objects[sphere_id]->is_mirror) {
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
        Vector P_light = random_point;
        Vector N_light = (P_light - s.light->C).getNormalized();
        Vector direction_light = P_light - P;
        double d_light2 = direction_light.norm2();
        Ray ray_light(P + 0.01 * N, direction_light.getNormalized());
        Vector P_light_inter, N_light_inter;
        int sphere_id_light;
        double t_light;
        bool has_inter_light = s.intersect(ray_light, P_light_inter, N_light_inter, sphere_id_light, t_light);

        // Check if the light source is visible from the point of intersection
        if (has_inter_light && t_light * t_light < d_light2 && sphere_id_light == s.objects.size()-1) {
            // Calculate the illumination intensity
            Vector albedo = s.objects[sphere_id]->albedo;
            double cos_theta = std::max(0.0, dot(N, direction_light.getNormalized()));
            double cos_light_theta = std::max(0.0, dot(N_light, -direction_light.getNormalized()));
            double cos_light_phi = std::max(0.0, dot(N_light, -axe_OP.getNormalized()));
            double distance_light2 = d_light2;
            intensite_pix = (albedo / M_PI) * (s.intensity_light * cos_theta * cos_light_theta * cos_light_phi / distance_light2);
        }
        else {
            intensite_pix = Vector(0, 0, 0);
        }
        }
    }

    // Background color
    else {
        intensite_pix = Vector(0, 0, 0);
    }

    return intensite_pix;
}

int main() {
    int W = 800; // Image width
    int H = 800; // Image height
    int nb_samples = 200; // Number of samples per pixel
    int nb_bounces = 2; // Number of bounces for ray tracing

    Scene scene;
    // Add light source
    Sphere light(Vector(0, 0, 20), 1.5, Vector(1, 1, 1), false, false);
    scene.light = &light;
    scene.intensity_light = 10; // Intensity of the light source

    // Add objects to the scene
    Sphere sphere1(Vector(0, 0, 3), 1, Vector(0.9, 0.3, 0.3), false, false);
    Sphere sphere2(Vector(2, 0, 4), 1, Vector(0.3, 0.9, 0.3), true, false);
    Sphere sphere3(Vector(-2, 0, 4), 1, Vector(0.3, 0.3, 0.9), false, false);
    Triangle triangle1(Vector(-5, -2, 2), Vector(5, -2, 2), Vector(0, 5, 2), Vector(0.5, 0.5, 0.5), false, true);
    scene.addSphere(sphere1);
    scene.addSphere(sphere2);
    scene.addSphere(sphere3);
    scene.addTriangle(triangle1);

    // Create an image buffer
    std::vector<unsigned char> image(W * H * 3, 0);

    // Camera parameters
    Vector C(0, 0, 0); // Camera center
    double dist = 1; // Distance from the camera to the image plane
    Vector I(1, 0, 0); // Camera direction
    Vector J(0, 1, 0); // Up vector
    Vector K(0, 0, 1); // Right vector

    // Loop through each pixel
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color(0, 0, 0);

            // For each sample in the pixel
            for (int s = 0; s < nb_samples; s++) {
                // Calculate random offsets for antialiasing
                double u = (j + uniform(engine)) / W;
                double v = (i + uniform(engine)) / H;

                // Compute ray direction
                Vector D = I * dist + (u - 0.5) * K + (v - 0.5) * J;

                // Create the ray
                Ray r(C, D.getNormalized());

                // Accumulate color using ray tracing
                color += getColor(r, scene, nb_bounces);
            }

            // Average the colors
            color = color / nb_samples;

            // Convert color to integer values (0-255) and store in image buffer
            image[(i * W + j) * 3 + 0] = std::min(255.0, std::max(0.0, color[0] * 255));
            image[(i * W + j) * 3 + 1] = std::min(255.0, std::max(0.0, color[1] * 255));
            image[(i * W + j) * 3 + 2] = std::min(255.0, std::max(0.0, color[2] * 255));
        }
    }

    // Write the image buffer to a PNG file
    stbi_write_png("output.png", W, H, 3, &image[0], 0);

    return 0;
}
