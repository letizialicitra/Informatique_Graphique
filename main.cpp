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
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0]=x;
        coord[1]=y;
        coord[2]=z;
    }
    // Access operacoord[1]=y;coord[2]=z;
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
    Vector getNormalized() {
        Vector result = *this;
        result.normalize();
        return result;
    }

    Vector operator-() const {
        return Vector(-coord[0], -coord[1], -coord[2]);
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
    // Constructor to initialize the sphere with a center, radius, and surface color
    Sphere(const Vector& c, double r, const Vector& a, bool is_mirror = false, bool is_transp = false) : C(c), R(r), albedo(a), is_mirror(is_mirror), is_transparent(is_transp) {};
    Vector C;       // Center of the sphere
    double R;       // Radius of the sphere
    Vector albedo;  // Surface color of the sphere
    bool is_mirror;
    bool is_transparent;

    // Function to check for intersection between the sphere and a ray
bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const{
        // resout a*t*2 + b*t +c =c0

        double a = 1;
        double b = 2 * dot(d.u, d.O - C);
        double c = (d.O - C).norm2() - R * R;

        double delta = b * b - 4 * a * c;
        if (delta < 0) return false;
        double t1 = (-b - sqrt(delta)) / 2 * a;
        double t2 = (-b + sqrt(delta)) / 2 * a;

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


// Class representing a scene containing multiple spheres
class Scene {
public:
    // Default constructor
    Scene() {};

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
    std::vector<Sphere> objects;
    Vector position_light;
    double intensity_light;
};

 


Vector getColor(Ray &r, const Scene &s, int nbrebonds) {

    if (nbrebonds == 0) return Vector(0, 0, 0);

    Vector P, N;
    int sphere_id;
    double t;
    bool has_inter = s.intersect(r, P, N, sphere_id, t);

    Vector intensite_pix(0, 0, 0);
    if (has_inter) {

        if (s.objects[sphere_id].is_transparent) {
            double n1 = 1;
            double n2 = 1.3;
            Vector normale_pour_transparence(N);
            if (dot(r.u, N) > 0) { // on sort de la sphere
                n1 = 1.3;
                n2 = 1;
                normale_pour_transparence = - N;
            }
            
            double radical = 1 - sqr(n1 / n2) * (1 - sqr(dot(normale_pour_transparence, r.u)));
            if (radical > 0) {
                Vector direction_refraction = (n1 / n2) * (r.u - dot(r.u, normale_pour_transparence) * normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
                Ray rayon_refracte(P - 0.01 * normale_pour_transparence, direction_refraction);
                intensite_pix = getColor(rayon_refracte, s, nbrebonds - 1);
            }
        }
        else if (s.objects[sphere_id].is_mirror) {
            Vector direction_mirroir = r.u - 2 * dot(N, r.u) * N;
            Ray rayon_mirroir(P + 0.01 * N, direction_mirroir);
            intensite_pix = getColor(rayon_mirroir, s, nbrebonds - 1);

        }
        else {

            Ray ray_light(P + 0.01 * N, (s.position_light - P).getNormalized());
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersect(ray_light, P_light, N_light, sphere_id_light, t_light);
            double d_light2 = (s.position_light - P).norm2();
            if (has_inter_light && t_light * t_light < d_light2) {
                intensite_pix = Vector(0, 0, 0);
            }
            else {
                intensite_pix = s.objects[sphere_id].albedo * (s.intensity_light * std::max(0., dot((s.position_light - P).getNormalized(), N)) / (s.position_light - P).norm2());
            }
        }
    }
    return intensite_pix;
}


int main() {
    int W = 512;
    int H = 512;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));

   

    //with this 2 you can reflect the color of the top on the ground: look at the z parameter, they are inversed! And this works just having as y (or as x -> you will have left and right of the same color) 10 or -10
    //  s.addSphere(Sphere(Vector(0, -10, 55), 10., Vector(0, 0.971, 0),true));         // ball in the center
    // s.addSphere(Sphere(Vector(0, 10, -55), 10., Vector(0, 0.971, 0),true));         // ball in the center
     // Create the scene and add spheres to it
    
   
    Scene s;
    Sphere s1(Vector(-15, 0, -55), 10., Vector(1,0,0),false,true);
    Sphere s1bis(Vector(15, 0, -55), 10., Vector(1,0,0),true);
    Sphere s2(Vector(0, -1000, 0), 960., Vector(0.0, 0.4, 0.14));   // from below
    Sphere s3(Vector(0, 1000, 0), 960., Vector(0.2, 0.2, 0.9));    // from up
    Sphere s4(Vector(-1000, 0, 0), 965., Vector(0.0, 0.2, 0.9));   // from left
    Sphere s5(Vector(1000, 0, 0), 965., Vector(0.9, 0.5, 0.7));    // from right
    Sphere s6(Vector(0, 0, -1000), 900., Vector(0.2, 0.1, 0.1));   // background !it is important to have the radius no more than 900, otherwise it is not possible to see the transparency
    // for example, if you put s6 with radius equal to 940, the spheres at the center will be too close to the background, so the shadow will cover the trasnparrency effect

 
    s.addSphere(s1);
    s.addSphere(s1bis);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    s.position_light = Vector(0, 30, 20);
    s.intensity_light = 1E9;
    Vector camera(0, 0 , 55);    // Set up the camera position
  
   
    // Create an image buffer
    std::vector<unsigned char> image(W * H * 3, 0);

    // int objectId;
    // double best_t;

    // Ray tracing loop
#pragma omp parallel for 
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // Calculate the ray direction based on the camera and image coordinates
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Ray r(camera, u);

            Vector color = getColor(r,s,5);
           
            // Apply gamma correction and store the color values in the image buffer
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2)));  // BLUE
        }
    }

    // Save the rendered image to a file
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
