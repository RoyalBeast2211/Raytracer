#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#if defined __linux__ || defined __APPLE__
#else
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

// A generic 3D vector class templated on numeric type T
template <typename T>
class Vec3
{
public:
    // Components of the vector
    T x, y, z;

    // Default constructor: initializes to zero vector
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}

    // Constructor: all components set to the same value
    Vec3(T xx) : x(xx), y(xx), z(xx) {}

    // Constructor: initializes each component separately
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    // Normalize the vector to length 1
    Vec3 &normalize()
    {
        // Compute the squared length of the vector
        T nor2 = length2();
        // Avoid division by zero
        if (nor2 > 0)
        {
            // Compute the inverse of the length
            T invNor = 1 / sqrt(nor2);
            // Scale each component
            x *= invNor;
            y *= invNor;
            z *= invNor;
        }
        // Return reference to this vector
        return *this;
    }

    // Scalar multiplication (vector * scalar)
    Vec3<T> operator*(const T &f) const
    {
        return Vec3<T>(x * f, y * f, z * f);
    }

    // Component-wise multiplication (vector * vector)
    Vec3<T> operator*(const Vec3<T> &v) const
    {
        return Vec3<T>(x * v.x, y * v.y, z * v.z);
    }

    // Dot product (scalar result)
    T dot(const Vec3<T> &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    // Vector subtraction
    Vec3<T> operator-(const Vec3<T> &v) const
    {
        return Vec3<T>(x - v.x, y - v.y, z - v.z);
    }

    // Vector addition
    Vec3<T> operator+(const Vec3<T> &v) const
    {
        return Vec3<T>(x + v.x, y + v.y, z + v.z);
    }

    // Compound addition: adds another vector to this one in place
    Vec3<T> &operator+=(const Vec3<T> &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    // Compound component-wise multiplication in place
    Vec3<T> &operator*=(const Vec3<T> &v)
    {
        x *= v.x;
        y *= v.y;
        z *= v.z;
        return *this;
    }

    // Unary minus: negates all components
    Vec3<T> operator-() const
    {
        return Vec3<T>(-x, -y, -z);
    }

    // Compute squared length of the vector (avoids sqrt)
    T length2() const
    {
        return x * x + y * y + z * z;
    }

    // Compute the length (magnitude) of the vector
    T length() const
    {
        return sqrt(length2());
    }

    // Friend function to allow easy printing of the vector
    friend std::ostream &operator<<(std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

// Type alias: Vec3f is Vec3<float>
typedef Vec3<float> Vec3f;

// Class representing a sphere in 3D space
class Sphere
{
public:
    // Center position of the sphere
    Vec3f center;

    // Sphere radius and radius squared (precomputed for efficiency)
    float radius, radius2;

    // Surface color (diffuse color) and emission color (light emitted by the sphere)
    Vec3f surfaceColor, emissionColor;

    // Transparency (0 = opaque) and reflectivity (0 = non-reflective)
    float transparency, reflection;

    // Constructor to initialize the sphere with position, radius, color, reflectivity, transparency, and emission
    Sphere(
        const Vec3f &c,          // center position
        const float &r,          // radius
        const Vec3f &sc,         // surface color
        const float &refl = 0,   // reflection coefficient (optional, default = 0)
        const float &transp = 0, // transparency coefficient (optional, default = 0)
        const Vec3f &ec = 0)     // emission color (optional, default = black)
        : center(c),
          radius(r),
          radius2(r * r), // precompute radius squared to avoid repeated multiplication
          surfaceColor(sc),
          emissionColor(ec),
          transparency(transp),
          reflection(refl)
    {
        // Empty body
    }

    // Compute a ray-sphere intersection using the geometric solution
    bool intersect(
        const Vec3f &rayorig, // origin of the ray
        const Vec3f &raydir,  // direction of the ray (assumed normalized)
        float &t0,            // output: distance to first intersection point
        float &t1) const      // output: distance to second intersection point
    {
        // Vector from ray origin to sphere center
        Vec3f l = center - rayorig;

        // Project vector l onto the ray direction
        float tca = l.dot(raydir);

        // If projection is negative, sphere is behind the ray
        if (tca < 0)
            return false;

        // Compute squared distance from sphere center to the ray
        float d2 = l.dot(l) - tca * tca;

        // If distance is larger than the radius, the ray misses the sphere
        if (d2 > radius2)
            return false;

        // Compute half-chord distance (thc)
        float thc = sqrt(radius2 - d2);

        // Compute intersection distances along the ray
        t0 = tca - thc;
        t1 = tca + thc;

        // There is an intersection
        return true;
    }
};

// This macro defines the maximum recursion depth for ray tracing
// (i.e., how many times a ray can bounce/reflection/refraction before stopping)
#define MAX_RAY_DEPTH 5

// Linear interpolation (mixing) function
// Parameters:
//   a   - first value
//   b   - second value
//   mix - mixing factor in [0,1]
// Returns:
//   a value interpolated between a and b:
//     - if mix=0, returns a
//     - if mix=1, returns b
//     - in between gives weighted blend
float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

// This is the main ray tracing function.
// It computes the color returned by a ray shot from 'rayorig' in direction 'raydir'.
// It recursively handles reflection, refraction, and shading.
Vec3f trace(
    const Vec3f &rayorig,               // Origin of the ray
    const Vec3f &raydir,                // Direction of the ray (should be normalized)
    const std::vector<Sphere> &spheres, // All spheres in the scene
    const int &depth)                   // Current recursion depth
{
    // Optional sanity check to ensure direction is normalized
    // if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;

    // Distance to the closest intersection so far
    float tnear = INFINITY;

    // Pointer to the sphere that is closest and hit by the ray
    const Sphere *sphere = NULL;

    // Find intersection of this ray with all spheres in the scene
    for (unsigned i = 0; i < spheres.size(); ++i)
    {
        float t0 = INFINITY, t1 = INFINITY;
        // Check if the ray intersects sphere i
        if (spheres[i].intersect(rayorig, raydir, t0, t1))
        {
            // If the nearest intersection is behind the origin, use the farther one
            if (t0 < 0)
                t0 = t1;
            // Update nearest intersection if this is the closest so far
            if (t0 < tnear)
            {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }

    // If no intersection was found, return background color (gray here)
    if (!sphere)
        return Vec3f(2);

    // Color accumulator for this ray
    Vec3f surfaceColor = 0;

    // Compute the intersection point
    Vec3f phit = rayorig + raydir * tnear;

    // Compute the normal at the intersection point
    Vec3f nhit = phit - sphere->center;
    nhit.normalize();

    // Bias to avoid self-intersection ("shadow acne")
    float bias = 1e-4;

    // Flag to track whether we are inside the sphere
    bool inside = false;

    // Flip the normal if we are inside the sphere
    if (raydir.dot(nhit) > 0)
    {
        nhit = -nhit;
        inside = true;
    }

    // Handle reflective and/or transparent materials recursively
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH)
    {
        // Compute facing ratio: how much the ray is aligned with the normal
        float facingratio = -raydir.dot(nhit);

        // Compute Fresnel effect: more reflection when looking at shallow angles
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);

        // Compute the reflection ray direction
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        refldir.normalize();

        // Trace reflection ray recursively
        Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);

        Vec3f refraction = 0;

        // If the sphere is transparent, compute refraction
        if (sphere->transparency)
        {
            float ior = 1.1;                      // Index of refraction
            float eta = (inside) ? ior : 1 / ior; // Inside or outside the sphere
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);

            // If k >=0, compute the refraction direction
            Vec3f refrdir = raydir * eta + nhit * (eta * cosi - sqrt(k));
            refrdir.normalize();

            // Trace refraction ray recursively
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }

        // Combine reflection and refraction contributions
        surfaceColor = (reflection * fresneleffect +
                        refraction * (1 - fresneleffect) * sphere->transparency) *
                       sphere->surfaceColor;
    }
    else
    {
        // Diffuse shading: for each light source, add contribution if not in shadow
        for (unsigned i = 0; i < spheres.size(); ++i)
        {
            // Only consider spheres with emission (light sources)
            if (spheres[i].emissionColor.x > 0)
            {
                Vec3f transmission = 1;
                Vec3f lightDirection = spheres[i].center - phit;
                lightDirection.normalize();

                // Check for shadows: cast a shadow ray to the light
                for (unsigned j = 0; j < spheres.size(); ++j)
                {
                    if (i != j)
                    {
                        float t0, t1;
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1))
                        {
                            // There is an occlusion: point is in shadow
                            transmission = 0;
                            break;
                        }
                    }
                }

                // Compute Lambertian shading (cosine of angle between normal and light)
                surfaceColor += sphere->surfaceColor * transmission *
                                std::max(float(0), nhit.dot(lightDirection)) *
                                spheres[i].emissionColor;
            }
        }
    }

    // Return the final color: shading + emission (if any)
    return surfaceColor + sphere->emissionColor;
}

// Main rendering function.
// This function generates an image by tracing a primary ray through each pixel.
// For each pixel:
//   - Compute the camera ray
//   - Trace the ray
//   - Store the resulting color
// Finally, it saves the image in PPM format.
void render(const std::vector<Sphere> &spheres)
{
    // Image resolution
    unsigned width = 800, height = 800;

    // Allocate memory to store pixel colors
    Vec3f *image = new Vec3f[width * height];
    Vec3f *pixel = image; // pointer to the current pixel

    // Precompute inverses for efficiency
    float invWidth = 1 / float(width);
    float invHeight = 1 / float(height);

    // Field of view in degrees
    float fov = 80;

    // Aspect ratio (width / height)
    float aspectratio = width / float(height);

    // Compute camera angle (convert FOV to radians and apply tangent)
    float angle = tan(M_PI * 0.5 * fov / 180.);

    // Loop over image rows (top to bottom)
    for (unsigned y = 0; y < height; ++y)
    {
        // Loop over columns (left to right)
        for (unsigned x = 0; x < width; ++x, ++pixel)
        {
            // Compute normalized device coordinates for pixel center
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;

            // Create ray direction from camera origin through pixel
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();

            // Trace the ray and store the resulting color
            *pixel = trace(Vec3f(0), raydir, spheres, 0);
        }
    }

    // Save the image to a PPM file
    // (binary P6 format, RGB bytes, 0-255)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);

    // Write PPM header
    ofs << "P6\n"
        << width << " " << height << "\n255\n";

    // Write pixel data
    for (unsigned i = 0; i < width * height; ++i)
    {
        // Clamp color components to [0,1], scale to [0,255]
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255)
            << (unsigned char)(std::min(float(1), image[i].y) * 255)
            << (unsigned char)(std::min(float(1), image[i].z) * 255);
    }

    ofs.close(); // close file

    // Clean up memory
    delete[] image;
}

// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene by calling the render() function.
int main(int argc, char **argv)
{
    // Seed the random number generator (not really used in this example)
    srand48(13);

    // Create a list of spheres to define the scene
    std::vector<Sphere> spheres;

    // Add spheres to the scene
    // Each sphere is defined by:
    //   - center position
    //   - radius
    //   - surface color
    //   - reflectivity
    //   - transparency
    //   - emission color (used for lights)

    // Large sphere to serve as the ground plane
    spheres.push_back(Sphere(
        Vec3f(0.0, -10004, -20), // Center far below the scene to act as a floor
        10000,                   // Huge radius
        Vec3f(0.20, 0.20, 0.20), // Grayish color
        0,                       // No reflectivity
        0.0));                   // Opaque

    // Main central sphere (colored, reflective and partly transparent)
    spheres.push_back(Sphere(
        Vec3f(0.0, 0, -20),      // Centered in front of camera
        4,                       // Radius
        Vec3f(1.00, 0.32, 0.36), // Redish color
        1,                       // Reflective
        0.5));                   // Semi-transparent

    // Small sphere on the right
    spheres.push_back(Sphere(
        Vec3f(5.0, -1, -15),     // Position
        2,                       // Radius
        Vec3f(0.90, 0.76, 0.46), // Yellowish color
        1,                       // Reflective
        0.0));                   // Opaque

    // Large sphere further back
    spheres.push_back(Sphere(
        Vec3f(5.0, 0, -25),      // Position
        3,                       // Radius
        Vec3f(0.65, 0.77, 0.97), // Blueish color
        1,                       // Reflective
        0.0));                   // Opaque

    // Sphere on the left
    spheres.push_back(Sphere(
        Vec3f(-5.5, 0, -15),     // Position
        3,                       // Radius
        Vec3f(0.90, 0.90, 0.90), // White color
        1,                       // Reflective
        0.0));                   // Opaque

    // Light source (emissive sphere)
    spheres.push_back(Sphere(
        Vec3f(0.0, 20, -30),     // Position high above
        3,                       // Radius
        Vec3f(0.00, 0.00, 0.00), // No diffuse color (black)
        0,                       // Non-reflective
        0.0,                     // Non-transparent
        Vec3f(3)));              // Emission color (bright white light)

    // Render the scene to an image file
    render(spheres);

    // Exit the program
    return 0;
}

// tinker with values of spheres to generate different image