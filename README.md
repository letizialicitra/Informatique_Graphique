# Informatique_Graphique
The aim of the course is to implement path tracing: it is a rendering technique used to generate realistic images by simulating the behavior of light within a three-dimensional sceneDifferent functionalities will be implemented at each BE.

# BE-0 
Before making any modifications, this code is the starting point.

```cpp
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static inline double sqr(double x) { return x * x; }

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    Vector& operator+=(const Vector& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

int main() {
    int W = 512;
    int H = 512;

    Vector center(0.2, 0.1, 0.);

    std::vector<unsigned char> image(W*H * 3, 0);
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            Vector v(j / (double)W - 0.5, i / (double)H - 0.5, 0.);
            double gaussianVal = exp(-(v-center).norm2()/(2*sqr(0.2)));

            image[(i*W + j) * 3 + 0] = 127* gaussianVal;   // RED
            image[(i*W + j) * 3 + 1] = 50 * gaussianVal;  // GREEN
            image[(i*W + j) * 3 + 2] = 255 * gaussianVal;  // BLUE
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
```

- `Vector` class represents a 3D vector and provides methods for vector operations such as addition, subtraction, multiplication, and calculating the squared norm.
- The `main()` function initializes the width `W` and height `H` of the image and creates a vector to store the image data.
- The nested loops iterate over each pixel of the image.
- Inside the loops, it calculates the Gaussian value at each pixel using the formula `exp(-(v-center).norm2()/(2*sqr(0.2)))`, where `v` is the vector representing the current pixel position, `center` is the center of the Gaussian distribution, and `0.2` is the standard deviation.
- The Gaussian value is then scaled and assigned to the red, green, and blue components of the pixel.
- Finally, the image is saved to a PNG file using the `stbi_write_png` function.

This code generates an image with a Gaussian distribution centered at `(0.2, 0.1, 0)`, where the intensity of each pixel is determined by the Gaussian function. The resulting image will have a Gaussian-like pattern, with brighter regions near the center and fading towards the edges.

This is the image produced:
![image_0](image_0.png)

# BE-1 12/01/24

## Section 1.0
First of all, we will introduce some objects, such as the Classes Ray and Sphere.
```cpp
// Class Ray -> // Class representing a ray in 3D space
class Ray {
public:
    Vector O; // Origin of the ray
    Vector u; // Direction of the ray

    // Constructor to initialize the ray with an origin and direction
    Ray(const Vector& O, const Vector& u) : O(O), u(u) {}
};
```
``` cpp
// Class representing a sphere in 3D space
class Sphere {
public:
    // Constructor to initialize the sphere with a center, radius, and surface color
    Sphere(const Vector& c, double r) : C(c), R(r) {};
    Vector C;       // Center of the sphere
    double R;       // Radius of the sphere
};
 ```

 After that, we have to create the image in the main and see if for every pixel of the image, the rayons from the camera will have intersections with the sphere.

 ```cpp
 int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));
    Vector camera(0, 0 , 55);    // Set up the camera position
    std::vector<unsigned char> image(W * H * 3, 0);  // Create an image buffer

    Sphere s1(Vector(0, 0, -55), 20.);  // Define a sphere
    
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);  // Calculate the ray direction based on the camera and image coordinates
            u.normalize();
            Ray r(Vector(0,0,0), u);
            bool intersection = s1.intersect(r);
            image[(i*W + j) * 3 + 0] = intersection? 255:0;   // RED
            image[(i*W + j) * 3 + 1] = intersection? 255:0;  // GREEN
            image[(i*W + j) * 3 + 2] = intersection? 255:0;  // BLUE
        }
    }
    stbi_write_png("image_1_0.png", W, H, 3, &image[0], 0);  // Save the image to a file

    return 0;
}
 ```
Now, we have to create the function 'intersect' for the Class Sphere. According to the result of delta, we can have or not have intersections.
 ```cpp
// Function to check for intersection between the sphere and a ray
bool intersect(const Ray& d) const{
    // Define coefficients of the quadratic equation for intersection
    double a = 1;
    double b = 2 * dot(d.u, d.O - C);
    double c = (d.O - C).norm2() - R * R;

    // Calculate the discriminant
    double delta = b * b - 4 * a * c;

    // If the discriminant is negative, no intersection
    if (delta < 0) 
        return false;

    // Calculate the solutions for t
    double t1 = (-b - sqrt(delta)) / (2 * a);
    double t2 = (-b + sqrt(delta)) / (2 * a);

    if (t2 > 0) 
        return true;

    return false;
}
 ```

 At the end of this section, we will obtain a white sphere on a black background. The color white is given cause of the rayons from the camera that have intersection with the area defined by the sphere. 

 ![image_1_0](image_1_0.png)

 ## Section 1.1
In this part, we will add the enlightenment Lambertien model.  We have to:
- add the source of light
- in case of intersection, take just the point closer to the light
- add the color (albedo) as feature of the Sphere.
  
This are the modification in the main.
```cpp
int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));
     Vector camera(0, 0 , 55);    // Set up the camera position
    // Create an image buffer
    std::vector<unsigned char> image(W * H * 3, 0);
    Sphere s1(Vector(0, 0, -55), 20., Vector(0.,1.,0.));

    Vector position_light = Vector(15, 70, -20);  // Set the position of the light source
    int intensity_light = 2E6;  // Set the intensity of the light source
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            // Calculate the ray direction based on the camera and image coordinates
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Vector P,N;
            Ray r(Vector(0,0,0), u);
            bool intersection =  s1.intersect(r,P,N);
            Vector intensity_pixel = Vector(0,0,0);
            if (intersection) {
                // Calculate the intensity of the pixel based on the Phong reflection model
                intensity_pixel = s1.albedo * intensity_light * std::max(0., dot(N, (position_light - P).getNormalized()))/(position_light - P).norm2();
                
            }
                // Clamp the pixel intensity values between 0 and 255
                image[(i*W + j) * 3 + 0] = std::min(255.,std::max(0., intensity_pixel[0]));   // RED
                image[(i*W + j) * 3 + 1] = std::min(255.,std::max(0., intensity_pixel[1]));  // GREEN
                image[(i*W + j) * 3 + 2] = std::min(255.,std::max(0., intensity_pixel[2]));  // BLUE
        }
    }
    stbi_write_png("image_1_1.png", W, H, 3, &image[0], 0);  // Save the rendered image

    return 0;
}

```

In the function intersect of the class Sphere, we had to add two Vectors, N (Normal) and P (position) in order to get the closer point of intersection. 
```cpp
bool intersect(const Ray& d, Vector& P, Vector& N) const{
        // resout a*t*2 + b*t +c =c0
        double t;
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

```
This is the image created: 
![image_1_1](image_1_1.png)

## Section 1.2
In this section we'll deal with the intersection between the rayons and the Scene. The Scene is a class that contains a vector of Spheres. It contains also a function to see if there are intersection for each Sphere of the Scene with the rays of the camera. 

```cpp
// Class representing a scene containing multiple spheres
// Definition of the Scene class
class Scene {
public:
    // Default constructor
    Scene() {}; // Empty constructor, does nothing

    // Function to add a sphere to the scene
    void addSphere(const Sphere& sphere) {
        objects.push_back(sphere); // Adds the sphere to the list of objects in the scene
    }

    // Function to check for intersection between the scene and a ray
    bool intersect(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
        bool has_inter = false; // Flag to indicate if there has been an intersection
        min_t = 1E99; // Initial maximum value for the minimum intersection distance
        for (int i = 0; i < objects.size(); ++i) { // Iterates over all objects in the scene
            Vector localP, localN; // Temporary variables to store local point and normal
            double t; // Temporary variable to store the intersection distance
            bool local_has_inter = objects[i].intersect(d, localP, localN); // Checks intersection between the current object and the ray
            if (local_has_inter) { // If there is an intersection with the current object
                has_inter = true; // Set the intersection flag to true
                if (t < min_t) { // If the intersection distance is less than the minimum distance found so far
                    min_t = t; // Update the minimum distance
                    P = localP; // Store the intersection point
                    N = localN; // Store the normal at the intersection point
                    sphere_id = i; // Store the ID of the object
                }
            }
        }
        return has_inter; // Returns true if there has been an intersection, otherwise false
    }

    std::vector<Sphere> objects; // Vector containing the objects present in the scene
};

```
We have to modify the main as well, adding a Scene and different spheres at it, to build the walls of our Scene.

```cpp
int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));
     Vector camera(0, 0 , 55);    // Set up the camera position
    // Create an image buffer
    std::vector<unsigned char> image(W * H * 3, 0);

     Scene s;
    
    Sphere s1(Vector(0, 0, -55), 12., Vector(0, 0.8, 0)); // center
    Sphere s2(Vector(0, -1000, 0), 985., Vector(1, 0.4, 0.14)); // down
    Sphere s3(Vector(0, 1000, 0), 910., Vector(0, 0, 0.8)); // up
    Sphere s4(Vector(-1000, 0, 0), 920., Vector(0.0, 0.5, 0.3)); // left
    Sphere s5(Vector(1000 - 15, 0, 0), 920., Vector(1, 1.0, 0)); // right
    Sphere s6(Vector(0, 0, -1000), 940., Vector(0.5, 0, 0)); // back

    s.addSphere(s1);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    Vector position_light = Vector(10,60,-20);
    int intensity_light = 4E6;
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            // Calculate the ray direction based on the camera and image coordinates
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Vector P,N;
            int sphere_id;
            Ray r(Vector(0,0,0), u);
            bool intersection =  s.intersect(r,P,N,sphere_id);
            Vector intensity_pixel = Vector(0,0,0);
            if (intersection) {
                intensity_pixel = s.objects[sphere_id].albedo * intensity_light * std::max(0., dot(N, (position_light - P).getNormalized()))/(position_light - P).norm2();
                
            }
                image[(i*W + j) * 3 + 0] = std::min(255.,std::max(0., intensity_pixel[0]));   // RED
                image[(i*W + j) * 3 + 1] = std::min(255.,std::max(0., intensity_pixel[1]));  // GREEN
                image[(i*W + j) * 3 + 2] = std::min(255.,std::max(0., intensity_pixel[2]));  // BLUE
        }
    }
    stbi_write_png("image_1_2.png", W, H, 3, &image[0], 0);

    return 0;
}

```

This is the image obtained at the end of the section:
![image_1_2](image_1_2.png)
# BE-2 19/01/24
## Section2.1 
In this part, we will introduce the shadowing. In particular, we will introduce, in case of intersection between the source of light and the Scene, another ray. In case this ray will have intersection with the Sphere and the distance between the sphere and the groundfloor is smaller than the distance between the the source of light and the groundfloor, we will put that pixel black.  

More over, we can add the gamma correction in order to have a major precision about the intensity of the color. For that, we use the function std::pow(intensity_color[0],1/2.2);
```cpp
int main() {
    int W = 512;
    int H = 512;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));
    Scene s;
    s.addSphere(Sphere(Vector(0, 0, 0), 10., Vector(0, 0.971, 0)));         // ball in the center
    s.addSphere(Sphere(Vector(-1000, 0, 0), 955., Vector(0.0, 0.2, 0.9)));   // from left
    s.addSphere(Sphere(Vector(1000, 0, 0), 965., Vector(0.9, 0.5, 0.7)));    // from right
    s.addSphere(Sphere(Vector(0, -1000, 0), 990., Vector(0.2, 0.4, 0.14)));   // from below
    s.addSphere(Sphere(Vector(0, 1000, 0), 950., Vector(0.2, 0.2, 0.1)));    // from up
    s.addSphere(Sphere(Vector(0, 0, -1000), 940., Vector(0.2, 0.1, 0.1)));   // background
    s.addSphere(Sphere(Vector(0, 0, 1000), 940., Vector(0.3, 0.4, 0.1)));    // background

    Vector camera(0, 0 , 55);    // Set up the camera position
    Vector light(-15, 30, 40);   // Set up the light position
    double intensity = 1E9;

    std::vector<unsigned char> image(W * H * 3, 0);

    int objectId;
    double best_t;

#pragma omp parallel for private(objectId, best_t)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Ray r(camera,u);
            Vector P,N;
            int sphere_id;
            double t;
            bool has_intersection = s.intersect(r, P, N, sphere_id, t);
            Vector color(0, 0, 0);

            if (has_intersection) {
                Ray ray_light(P + 0.01 * N, (light - P).getNormalized());
                Vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_intersection_light = s.intersect(ray_light, P_light, N_light, sphere_id_light, t_light);
                double d_light_squared = (light - P).norm2();

                if (has_intersection_light && t_light * t_light < d_light_squared) {
                    color = Vector(0, 0, 0); // Shadow color -> 000 is black, indicating shadow
                } else {
                    color = s.objects[sphere_id].albedo * (intensity * std::max(0., dot((light - P).getNormalized(), N)) / d_light_squared);
                }
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));   // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));   // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));   // BLUE
        }
    }

    stbi_write_png("image_2_1.png", W, H, 3, &image[0], 0);

    return 0;
}

```
This is the result obtained: it is possible to notice the shadow below the ball in the center, that goes onto the floor.qa
![image_2_1](image_2_1.png)

## Section2.2 Surface mirror

In this section we aim to introduce the concept of mirror surface. What we have to do is:
- create a function tha compute recursively the color to apply, called getColor();
- the number of calls of the function is determined by a parameter called number_of_rebonds;
- we called this function in the main; it is also a way to have a better style of coding; 
- we add a parameter inside the Sphere Class to know if that sphere is mirror or not.


```cpp
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


```

In the main:
```cpp
 for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Ray r(camera,u);
            Vector color = getColor(r,s,5);
            // Apply gamma correction and store the color values in the image buffer
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));   // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));   // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));   // BLUE
        }
    }

```
![image_2_2](image_2_2.png)


## Section2.2 Surface transparent
In this section we'll see hot to implement a trasnparent surface. The modifications that we have to do are the following:
- add a boolean as parameter of the Sphere to know if it is transparent or not;
- modify the function getColor, in case we have the boolean set to 1, we add to compute the rifraction of the ray and compute the tangent and the normal component
- add another sphere at the center of the image in order to see the difference between trasnaprency and mirroring.

```cpp
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

```

![image_2_3](image_2_3.png)
# BE-3 02/02/24

# BE-4 09/02/24

```cpp


```