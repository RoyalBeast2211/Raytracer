# Simple C++ Raytracer

This project is a minimalist raytracer written in a single C++ file. It's designed as an educational tool to understand the fundamental principles of computer graphics, rendering, and C++ programming. The code implements core ray tracing concepts to generate an image of a scene composed of spheres.

## Features

- **3D Vector Math:** A `Vec3` class for vector and color operations.
- **Ray-Sphere Intersection:** Geometric solution for finding intersections between rays and spheres.
- **Recursive Ray Tracing:** To handle reflections and refractions up to a maximum depth.
- **Shading and Lighting:**
  - Diffuse (Lambertian) shading.
  - Hard shadows.
  - Emissive materials for light sources.
- **Materials:**
  - Configurable surface color.
  - Reflection (mirrors).
  - Refraction (glass/transparency) with Fresnel effect.
- **Image Output:** Renders the scene to a `.ppm` image file.

## Prerequisites

To compile and run this project, you need a C++ compiler that supports C++11 or later. `g++` is recommended on Linux.

- A C++ compiler (e.g., GCC, Clang, MSVC)

## How to Compile and Run

1.  **Clone the repository:**

    ```bash
    git clone <your-repository-url>
    cd <repository-directory>
    ```

2.  **Compile the code:**
    Open your terminal and use the following command to compile `raytracer.cpp`. The `-O3` flag enables optimizations, which is highly recommended for a raytracer.

    ```bash
    g++ -o raytracer raytracer.cpp -O3 -std=c++11
    ```

3.  **Run the executable:**
    This will start the rendering process.

    ```bash
    ./raytracer
    ```

4.  **View the Output:**
    The program will generate an image file named `untitled.ppm` in the same directory. You can view this file using an image viewer that supports the PPM format, such as GIMP, IrfanView, or Preview on macOS.

## Customizing the Scene

You can easily modify the scene by changing the sphere properties inside the `main` function in `raytracer.cpp`.

The `Sphere` constructor is defined as:

```cpp
Sphere(position, radius, surface_color, reflectivity, transparency, emission_color)
```

### Example: Changing a Sphere's Color

To change the central sphere's color to blue, find its definition in `main()` and modify the `surfaceColor` vector.

**Original:**

```cpp

    // Main central sphere (colored, reflective and partly transparent)
    spheres.push_back(Sphere(
        Vec3f(0.0, 0, -20),      // Centered in front of camera
        4,                       // Radius
        Vec3f(1.00, 0.32, 0.36), // Redish color
        1,                       // Reflective
        0.5));                   // Semi-transparent
```

**Modified (to blue):**

```cpp

    // Main central sphere (colored, reflective and partly transparent)
    spheres.push_back(Sphere(
        Vec3f(0.0, 0, -20),      // Centered in front of camera
        4,                       // Radius
        Vec3f(0.2, 0.2, 1.0),    // Blue color
        1,                       // Reflective
        0.5));                   // Semi-transparent

```

Feel free to add, remove, or modify the spheres to create your own unique scenes!

## Code Structure

The entire logic is contained within `raytracer.cpp`:

- **`Vec3<T>`:** A template class for 3D vectors, used for points, directions, and colors.
- **`Sphere`:** A class to define a sphere with its physical and material properties.
- **`trace()`:** The core recursive function that traces a single ray and calculates the color it sees by handling intersections, reflections, and refractions.
- **`render()`:** Sets up the camera, iterates through each pixel of the image, casts a primary ray using `trace()`, and saves the final image to a PPM file.
- **`main()`:** The entry point of the program. This is where the scene (a vector of `Sphere` objects) is defined and the `render()` function is called.
