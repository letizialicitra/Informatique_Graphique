#include <vector>
#include <cmath>
#include <algorithm>

// Include STB image write and read implementations
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Define PI for later use
#define M_PI 3.14159265358979323846

// Function to calculate the square of a number
static inline double sqr(double x) { return x * x; }

// Vector class to represent 3D vectors
class Vector {
public:
    // Constructor with default values
    explicit Vector(double x = 0, double y = 0, double z = 0) : coord{x, y, z} {}

    // Access operators for vector components
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    // Overloaded addition assignment operator
    Vector& operator+=(const Vector& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    // Calculate the squared norm of the vector
    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    // Normalize the vector
    void normalize() {
        double n = sqrt(norm2());
        coord[0] /= n;
        coord[1] /= n;
        coord[2] /= n;
    }

    // Get the normalized vector without modifying the original vector
    Vector getNormalized() const {
        Vector result = *this;
        result.normalize();
        return result;
    }

    // Array to store the vector components
    double coord[3];
};

// Overloaded addition operator for vectors
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

// Overloaded subtraction operator for vectors
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// Overloaded multiplication operator for vector and scalar
Vector operator*(const Vector& a, double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

// Overloaded multiplication operator for scalar and vector
Vector operator*(double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

// Dot product of two vectors
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Class representing a ray in 3D space
class Ray {
public:
    Vector O; // Origin of the ray
    Vector u; // Direction of the ray

    // Constructor to initialize the ray with an origin and direction
    Ray(const Vector& O, const Vector& u) : O(O), u(u) {}
};

// Class representing a sphere in 3D space
class Sphere {
public:
    Vector C;       // Center of the sphere
    double R;       // Radius of the sphere
    Vector albedo;  // Surface color of the sphere
    bool is_mirror;
    // Function to check for intersection between the sphere and a ray
    bool intersect(const Ray& r, Vector& P, Vector& N, double& tlocal) const {
        // Coefficients of the quadratic equation for intersection
        double a = 1;
        double b = 2 * dot(r.u, r.O - C);
        double c = (r.O - C).norm2() - sqr(R);

        // Discriminant of the quadratic equation
        double delta = b * b - 4 * a * c;

        // If the discriminant is negative, there is no intersection
        if (delta < 0) return false;

        // Calculate the square root of the discriminant
        double sqrtdelta = sqrt(delta);

        // Calculate the two possible intersection points
        double t1 = (-b - sqrtdelta) / (2 * a);
        double t2 = (-b + sqrtdelta) / (2 * a);

        // If both intersection points are behind the ray, there is no intersection
        if (t2 < 0) return false;

        // Choose the nearest intersection point
        double t = t1;
        if (t1 < 0) t = t2;

        // Calculate the intersection point and the normal vector at that point
        P = r.O + t * r.u;
        N = P - C;
        N.normalize();

        // Store the local parameter of the intersection
        tlocal = t;

        // Intersection found
        return true;
    }

    // Constructor to initialize the sphere with a center, radius, and surface color
    Sphere(const Vector& c, double r, const Vector& a, bool is_mirror = false) : C(c), R(r), albedo(a), is_mirror(is_mirror) {}
};

// Class representing a scene containing multiple spheres
class Scene {
public:
    // Default constructor
    Scene() {}

    // Function to add a sphere to the scene
    void addSphere(const Sphere& sphere) {
        objects.push_back(sphere);
    }

    // Function to check for intersection between the scene and a ray
    bool intersect(const Ray& r, Vector& P, Vector& N, int& objectId, double& best_t) const {
        bool has_intersection = false;
        best_t = 1E10; // Initialize with a large value

        // Loop through all the objects in the scene
        for (int i = 0; i < objects.size(); i++) {
            Vector Plocal, Nlocal;
            double tlocal;

            // Check for intersection with the current object
            if (objects[i].intersect(r, Plocal, Nlocal, tlocal)) {
                has_intersection = true;

                // Update the best intersection parameters if the current intersection is closer
                if (tlocal < best_t) {
                    best_t = tlocal;
                    objectId = i; // Store the index of the intersected object
                    P = Plocal;   // Store the intersection point
                    N = Nlocal;   // Store the normal vector at the intersection point
                }
            }
        }

        // Return true if there is any intersection in the scene
        return has_intersection;
    }

    // Vector to store the spheres in the scene
    std::vector<Sphere> objects;
    Vector position_light ;   // Set up the light position
    double intensity_light ;

};


Vector getColor(const Ray &r,const  Scene &s, int number_of_rebonds ) {

    if ( number_of_rebonds == 0 ) {
        return Vector(0,0,0);
    }
     Vector P, N;
            int sphere_id;
            double t;
            bool has_intersection = s.intersect(r, P, N, sphere_id, t);
            Vector color(0, 0, 0);


            if (has_intersection) {

              
                if(s.objects[sphere_id].is_mirror) {
                    Vector direction_mirror= r.u - 2*dot(N,r.u)*N;
                    Ray ray_mirror(P+0.001*N,direction_mirror);
                    color = getColor(ray_mirror, s, number_of_rebonds -1);
                } else {

                Ray ray_light(P + 0.01 * N, (s.position_light- P).getNormalized());
                Vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_intersection_light = s.intersect(ray_light, P_light, N_light, sphere_id_light, t_light);
                double d_light_squared = (s.position_light - P).norm2();

                if (has_intersection_light && t_light * t_light < d_light_squared) {
                    color = Vector(0, 0, 0); // Shadow color -> 000 è nero, fa 'ombra'
                } else {
                    color = s.objects[sphere_id].albedo * (s.intensity_light * std::max(0., dot((s.position_light - P).getNormalized(), N)) / d_light_squared);
                }
            }
            }
            return color;

}


int main() {
    int W = 512;
    int H = 512;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));

    // Create the scene and add spheres to it
    Scene s;

    //with this 2 you can reflect the color of the top on the ground: look at the z parameter, they are inversed! And this works just having as y (or as x -> you will have left and right of the same color) 10 or -10
    //  s.addSphere(Sphere(Vector(0, -10, 55), 10., Vector(0, 0.971, 0),true));         // ball in the center
    // s.addSphere(Sphere(Vector(0, 10, -55), 10., Vector(0, 0.971, 0),true));         // ball in the center
    s.addSphere(Sphere(Vector(-5, -15, -55), 10., Vector(0, 0.971, 0),true));         // ball in the center
    s.addSphere(Sphere(Vector(5, 15, -55), 10., Vector(0, 0.971, 0),true));         // ball in the center



    s.addSphere(Sphere(Vector(-1000, 0, 0), 965., Vector(0.0, 0.2, 0.9)));   // from left
    s.addSphere(Sphere(Vector(1000, 0, 0), 965., Vector(0.9, 0.5, 0.7)));    // from right
    s.addSphere(Sphere(Vector(0, -1000, 0), 960., Vector(0.0, 0.4, 0.14)));   // from below
    s.addSphere(Sphere(Vector(0, 1000, 0), 960., Vector(0.2, 0.2, 0.9)));    // from up
    s.addSphere(Sphere(Vector(0, 0, -1000), 940., Vector(0.2, 0.1, 0.1)));   // background
    s.addSphere(Sphere(Vector(0, 0, 1000), 940., Vector(0.3, 0.4, 0.1)));    // background
    s.position_light = Vector(-15, 30, 40);   // Set up the light position
    s.intensity_light = 1E9;
    Vector camera(0, 0 , 55);    // Set up the camera position
  
   
    // Create an image buffer
    std::vector<unsigned char> image(W * H * 3, 0);

    int objectId;
    double best_t;

    // Ray tracing loop
#pragma omp parallel for private(objectId, best_t)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // Calculate the ray direction based on the camera and image coordinates
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Ray r(camera, u);

            Vector color = getColor(r,s,5);
           
            // Apply gamma correction and store the color values in the image buffer
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));   // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));   // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));   // BLUE
        }
    }

    // Save the rendered image to a file
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
