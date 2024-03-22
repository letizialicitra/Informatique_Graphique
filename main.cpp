#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include "library.h"

Geometry::Geometry(const char * obj, double scaling, const Vector & offset, const Vector & color, bool mirror, bool transparent) {};

int main() {
    int W = 512;
    int H = 512;
    const int number_of_rays = 80;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));

    Scene scene;

    Sphere sphere_light(Vector(-5, 30, 40), 15, Vector(1.,1.,1.));

    Sphere s1(Vector(0, 0, -55), 10., Vector(1, 0, 0));
    // Sphere s1b(Vector(-15, 0, -35), 10., Vector(1, 1, 0.2),false,true);
    // Sphere s1c(Vector(15, 0, -75), 10., Vector(1, 0, 1),true);
    
    Sphere sphere1(Vector(0, 0, 3), 1., Vector(0.9, 0.3, 0.3), false, false);
    Sphere sphere2(Vector(2, 0, 4), 1., Vector(0.3, 0.9, 0.3), true, false);
    Sphere sphere3(Vector(-2, 0, 4), 1., Vector(0.3, 0.3, 0.9), false, false);
    Sphere s2(Vector(0, -1000, 0), 960., Vector(0.0, 0.4, 0.14));
    Sphere s3(Vector(0, 1000, 0), 960., Vector(0.2, 0.2, 0.9));
    Sphere s4(Vector(-1000, 0, 0), 965., Vector(0.0, 0.2, 0.9));
    Sphere s5(Vector(1000, 0, 0), 965., Vector(0.9, 0.5, 0.7));
    Sphere s6(Vector(0, 0, -1000), 900., Vector(0.2, 0.1, 0.1));

    //Geometry g1("cat.obj",100, Vector(0,0,-55), Vector(1.,1.,1.));
    Triangle t1(Vector(-5,-5,-20), Vector(5, -5, -20), Vector(0,5,-20), Vector(0.,1.,0.),false,false);
    // Sphere sphere1(Vector(0, 0, 3), 10., Vector(0.9, 0.3, 0.3), false, false);
    // Sphere sphere2(Vector(2, 0, 4), 10., Vector(0.3, 0.9, 0.3), true, false);
    // Sphere sphere3(Vector(-2, 0, 4), 10., Vector(0.3, 0.3, 0.9), false, false);
    //Triangle t1(Vector(-5, -2, 2), Vector(5, -2, 2), Vector(0, 5, 2), Vector(0., 1, 0.), false, false);
    scene.addSphere(sphere_light);

    scene.addSphere(s1);
    // scene.addSphere(s1b);
    // scene.addSphere(s1c);
    //scene.addGeometry(g1);
    scene.addSphere(sphere1); 
    scene.addSphere(sphere2);
    scene.addSphere(sphere3);
    scene.addSphere(s2);
    scene.addSphere(s3);
    scene.addSphere(s4);
    scene.addSphere(s5);
    scene.addSphere(s6);
    
    
    scene.addTriangle(t1);
    scene.light = &sphere_light;
    scene.intensity_light = 5E9;
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
                color += getColor(r, scene, 5) / number_of_rays;
            } 
            // Apply gamma correction and store the color values in the image buffer
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2))); // BLUE
        }
    }

    // Write the rendered image to file
    stbi_write_png("output.png", W, H, 3, &image[0], 0);

    return 0;
}
