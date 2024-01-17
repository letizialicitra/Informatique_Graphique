#include <vector>
#include <cmath>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265358979323846

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

    void normalize() {
        double n = sqrt(norm2());
        coord[0] /= n;
        coord[1] /= n;
        coord[2] /= n;
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
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray {
public:
    Vector O;
    Vector u;

    Ray(const Vector& O, const Vector& u) : O(O), u(u) {}
};

class Sphere {
public:
    Vector C;
    double R;
    Vector albedo;

    bool intersect(const Ray& r, Vector& P, Vector& N, double& tlocal) const {
        double a = 1;
        double b = 2 * dot(r.u, r.O - C);
        double c = (r.O - C).norm2() - sqr(R);
        double delta = b * b - 4 * a * c;

        if (delta < 0) return false;

        double sqrtdelta = sqrt(delta);
        double t1 = (-b - sqrtdelta) / (2 * a);
        double t2 = (-b + sqrtdelta) / (2 * a);

        if (t2 < 0) return false;

        double t = t1;
        if (t1 < 0) t = t2;

        P = r.O + t * r.u;
        N = P - C;
        N.normalize();

        tlocal = t;

        return true;
    }

    Sphere(const Vector& c, double r, const Vector& a) : C(c), R(r), albedo(a) {}
};

class Scene {
public:
    Scene() {}
    void addSphere(const Sphere& sphere) {
        objects.push_back(sphere);
    }

    bool intersect(const Ray& r, Vector& P, Vector& N, int& objectId) const {
        bool has_inter = false;
        double best_t = 1E10;
        for (int i = 0; i < objects.size(); i++) {
            Vector Plocal, Nlocal;
            double tlocal;
            if (objects[i].intersect(r, Plocal, Nlocal, tlocal)) {
                has_inter = true;
                if (tlocal < best_t) {
                    best_t = tlocal;
                    objectId = i;
                    P=Plocal;
                    N=Nlocal;
                }
            }
        }
        return has_inter;
    }


    std::vector<Sphere> objects;
};

int main() {
    int W = 512;
    int H = 512;

    double fov = 60 * M_PI / 180.0;
    double d = W / (2.0 * tan(fov / 2));
    Scene s;
    s.addSphere(Sphere(Vector(0, 0, 0), 10., Vector(0.3, 0.4, 0.1)));
    s.addSphere(Sphere(Vector(-1000, 0, 0), 940., Vector(0.3, 0.2, 0.9)));
    s.addSphere(Sphere(Vector(1000, 0, 0), 990., Vector(0.9, 0.5, 0.5)));
    s.addSphere(Sphere(Vector(0, -1000, 0), 940., Vector(0.8, 0.4, 0.9)));
    s.addSphere(Sphere(Vector(0, 1000, 0), 940., Vector(0.2, 0.2, 0.1)));
    s.addSphere(Sphere(Vector(0, 0, -1000), 940., Vector(0.3, 0.1, 0.9)));
    s.addSphere(Sphere(Vector(0, 0, 1000), 940., Vector(0.3, 0.4, 0.8)));

    Vector camera(0, 0, 55);
    Vector L(-10, 20, 40);
    double I = 1E10;

    std::vector<unsigned char> image(W * H * 3, 0);

    int objectId;
    double best_t;

#pragma omp parallel for private(objectId, best_t)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W / 2. + 0.5, -i + H / 2. - 0.5, -d);
            u.normalize();
            Ray r(camera, u);
            Vector P, N;
            if (s.intersect(r, P, N, objectId)) {
                Vector vecLum = L - P;
                double d2 = vecLum.norm2();
                vecLum.normalize();
                Vector color = s.objects[objectId].albedo * (I * std::max(0., dot(vecLum, N)) / (4. * M_PI * M_PI * d2));
                image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));   // RED
                image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));   // GREEN
                image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));   // BLUE
            }
        }
    }

    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
