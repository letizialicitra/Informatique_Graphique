#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cfloat>
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

    Vector albedo;
    bool is_mirror; 
    bool is_transparent; 

};

class BBox{

public:
    Vector pMin, pMax;
    
   
    BBox(const Vector &pMin, const Vector &pMax) : pMin(pMin), pMax(pMax) {}
    BBox() : pMin(Vector()), pMax(Vector()) {}
    void operator=(const BBox &b)
    {
        pMin = b.pMin;
        pMax = b.pMax;
    }

    Vector getMin() const
    {
        return pMin;
    }

    Vector getMax() const
    {
        return pMax;
    }

    bool intersect(const Ray &r) const {
        Vector invU = Vector(1 / r.u[0], 1 / r.u[1], 1 / r.u[2]);

        double tMinX = (pMin[0] - r.O[0]) * invU[0];
        double tMaxX = (pMax[0] - r.O[0]) * invU[0];
        if (tMinX > tMaxX)
            std::swap(tMinX, tMaxX);

        double tMinY = (pMin[1] - r.O[1]) * invU[1];
        double tMaxY = (pMax[1] - r.O[1]) * invU[1];
        if (tMinY > tMaxY)
            std::swap(tMinY, tMaxY);

        double tMinZ = (pMin[2] - r.O[2]) * invU[2];
        double tMaxZ = (pMax[2] - r.O[2]) * invU[2];
        if (tMinZ > tMaxZ)
            std::swap(tMinZ, tMaxZ);

        double t0 = std::max(tMinX, std::max(tMinY, tMinZ));
        double t1 = std::min(tMaxX, std::min(tMaxY, tMaxZ));

        return t0 < t1 && t1 > 0;
    }
};


class Geometry: public Object {
public:

    Geometry(const char* obj, double scaling, const Vector& offset, const Vector& color, bool mirror= false, bool transparent= false);
 
    std::vector<int> faceGroup;
    std::vector<int> faces;
    std::vector<int> normalIds;
    std::vector<int> uvIds;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    bool intersect(const Ray& d, Vector& P, Vector& N, double &t) const{
        {
        
        if(!bbox.intersect(d)) return false;
        
        bool has_inter = false;
        double min_t = 1E99;
        for (int i = 0; i < faces.size(); i+=3) {
            Vector A = vertices[faces[i]];
            Vector B = vertices[faces[i+1]];
            Vector C = vertices[faces[i+2]];
            Vector localP, localN;
            double local_t;
            bool local_has_inter = intersect(d, localP, localN, local_t);
            if (local_has_inter) {
                has_inter = true;
                if (local_t < min_t) {
                    min_t = local_t;
                    P = localP;
                    N = localN;
                }
            }
        }
        t = min_t;
        return has_inter;
    
    }
}

private:
    BBox bbox;
};



// Sphere class representing a sphere in 3D space
class Sphere :public Object {
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
    std::vector<const Object *> objects; // Collection of spheres
    Sphere *light;       //pointer to a spheric light source
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

    // Function to add a geometry to the scene
    void addGeometry(const Geometry& geometry) {
        objects.push_back(&geometry);
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
        
        intensite_pix = (s.intensity_light / (4*M_PI*d_light_squared) * std::max(0., dot(N, wi)) /dot(axe_OP,random_direction)) *s.objects[sphere_id]->albedo;
        }
        // Handling indirect illumination
     
        Vector random_dir = random_cos(N);
        Ray random_rayon(P + 0.001 * N, random_dir);
        intensite_pix += getColor(random_rayon, s, nbrebonds - 1) * s.objects[sphere_id]->albedo;
    }
}
    return intensite_pix;
}







class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};


class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vector &albedo, bool is_mirror = false, bool is_transparent = false) 
    {
        this->albedo = albedo;
        this->is_mirror = is_mirror;
        this->is_transparent = is_transparent;
    };



  

    void scale(double s)
    {
        center = Vector(0, 0, 0);
        for (size_t i = 0; i < vertices.size(); i++)
        {
            vertices[i] = vertices[i] * s;
            center = center + vertices[i];
        }

        if (vertices.size() != 0)
        {
            center = center * (1 / vertices.size());
        }
        bbox = BBox(bbox.getMin() * s, bbox.getMax() * s);
    }

    void translate(Vector t)
    {
        center = Vector(0, 0, 0);
        for (size_t i = 0; i < vertices.size(); i++)
        {
            vertices[i] += t;
            center = center + vertices[i];
        }

        if (vertices.size() != 0)
        {
            center = center * (1 / vertices.size());
        }
        bbox = BBox(bbox.getMin() + t, bbox.getMax() + t);
    }

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);

       
        center = Vector(0, 0, 0);
        double minx = INFINITY, miny = INFINITY, minz = INFINITY;
        double maxx = -INFINITY, maxy = -INFINITY, maxz = -INFINITY;
        for (size_t i = 0; i < vertices.size(); i++)
        {
            center = center + vertices[i];

            minx = std::min(minx, vertices[i][0]);
            miny = std::min(miny, vertices[i][1]);
            minz = std::min(minz, vertices[i][2]);

            maxx = std::max(maxx, vertices[i][0]);
            maxy = std::max(maxy, vertices[i][1]);
            maxz = std::max(maxz, vertices[i][2]);
        }
        center = center * (1 / vertices.size());
        bbox = BBox(Vector(minx, miny, minz), Vector(maxx, maxy, maxz));
    }

    bool intersect(const Ray &r, Vector &P, Vector &N, double &t) {
        if (!bbox.intersect(r))
            return false;

        bool has_inter = false;
        t = INFINITY;
        for (size_t i = 0; i < indices.size(); i++)
        {
            const Vector &A = vertices[indices[i].vtxi];
            const Vector &B = vertices[indices[i].vtxj];
            const Vector &C = vertices[indices[i].vtxk];

            Vector e1 = B - A;
            Vector e2 = C - A;
            Vector nLocal = cross(e1, e2);
            double uN = dot(r.u, nLocal);

            if (uN == 0)
                continue;

            double tLocal = dot(A - r.O, nLocal) / uN;
            if (tLocal < 0)
                continue;

            double b = dot(e2, cross(A - r.O, r.u)) / uN;
            if (b < 0 || b > 1)
                continue;
            double g = dot(e1, cross(r.u, A - r.O)) / uN;
            if (g < 0 || b + g > 1)
                continue;
            double a = 1 - b - g;

            if(b>=0 && g>=0 && b<=1 && g<=1 && a>=0 && tLocal>0){
                has_inter = true;
                if (tLocal < t)
                {
                    t = tLocal;
                    P = r.O + t * r.u;
                    N =nLocal;
                    N.normalize();
                }
            }
        }
        return has_inter;
    }

    Vector getCenter() const
    {
        return center;
    }

private:
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    Vector center;
    BBox bbox;
};






