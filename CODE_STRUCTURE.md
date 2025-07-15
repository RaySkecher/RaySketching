# Code Structure of `trace.c`

This document breaks down the structure and execution flow of the ray tracer implemented in `trace.c`.

## 1. Entry Point: `main()`

The program starts execution in the `main` function. Here's what it does:

1.  **`initialize_scene()`**: Sets up the Cornell Box scene, including walls, boxes, spheres, and their materials.
2.  **File Setup**: Opens `image.ppm` for writing and writes the PPM header (format, dimensions, max color value).
3.  **Pixel Iteration**: It iterates through each pixel of the output image (`width` x `height`), from top to bottom.
    -   For each pixel, it prints the "Scanlines remaining" to `stderr`.
4.  **Sampling Loop (Anti-aliasing)**: Inside the pixel loop, there's a loop that runs `SAMPLES` times for anti-aliasing.
    -   **Seed Generation**: A unique seed for the random number generator is created based on pixel coordinates and sample index.
    -   **UV Coordinates**: Calculates normalized `(u, v)` coordinates for the current sample within the pixel.
    -   **Camera Ray Generation**:
        -   `camera_const()`: Creates a camera for the scene.
        -   `camera_get_ray()`: Generates a ray from the camera through the pixel's `(u,v)` coordinates. This includes logic for depth of field if `aperture` is non-zero.
    -   **`color(ray)`**: This function is called to compute the color for the generated ray. This is the core of the ray tracing logic.
    -   **Accumulation**: The color returned from `color()` is added to a total color for the pixel (`tcol`).
5.  **Final Pixel Color Calculation**:
    -   **Averaging**: The accumulated color is divided by `SAMPLES` to get the average color.
    -   **Gamma Correction**: `sqrtf()` is applied to each color component (a simple form of gamma correction, gamma 2.0).
    -   **Clamping and Scaling**: The final color values are clamped to `[0, 255]` and converted to integers.
6.  **File Writing**: The final integer RGB values for the pixel are written to `image.ppm`.
7.  **Cleanup**: After the loops complete, it prints "Done." and closes the file.

## 2. Core Ray Tracing: `color(ray)`

This function calculates the color a given ray "sees" in the scene. It simulates light transport.

1.  **Initialization**: Initializes throughput color `col` to white `{1,1,1}` and `emitted` light to black `{0,0,0}`.
2.  **Bounce Loop**: It simulates the path of the ray for up to `MAX_BOUNCES`.
    -   **`world_hit(ray, ...)`**: Finds the closest object the ray intersects with.
        -   If no object is hit, the loop terminates, and it returns the accumulated `emitted` light.
    -   **Material Handling**:
        -   **Light Source**: If the ray hits a light source (`DIFFUSE_LIGHT`), it adds the light's color to the result. If it's the first bounce, it returns the light's color directly. Otherwise, it stops and returns what has been emitted so far. This is because we are doing explicit light sampling.
        -   **Lambertian Surface**:
            -   **`material_scatter()`**: Calculates a new scattered ray with a random direction from the hit point. It also returns the material's color (`attenuation`).
            -   The throughput `col` is updated by multiplying with the material's color.
            -   The input `ray` for the next bounce is replaced with the `scattered` ray.
    -   **Direct Light Sampling (Next Event Estimation)**:
        -   A random point on the light source (`box_light`) is chosen.
        -   A "shadow ray" is cast from the current hit point towards this point on the light.
        -   **`shadow_hit(...)`**: Checks if this shadow ray is blocked by any other object before it reaches the light.
        -   If the shadow ray is not blocked, the light's contribution is calculated (considering distance, angles, and area of the light) and added to the `emitted` color, weighted by the current `col` throughput.
3.  **Return Value**: The function returns the accumulated `emitted` light from direct light sampling over all bounces.

## 3. Scene and Intersections

### `world_hit()`

This function orchestrates intersection tests for a ray against all objects in the scene.

-   It calls intersection functions for each object (`box_intersect`, `sphere_hit`).
-   It keeps track of the closest intersection (`closest_so_far`).
-   If an intersection is found, it populates a `hit_record` struct with information about the hit (point, normal, material) and returns `true`.

### `shadow_hit()`

A specialized version of `world_hit` used for shadow rays.

-   It only needs to know *if* there is an intersection, not which one is closest or detailed hit info.
-   It calls specialized intersection functions (`box_intersect_shadow`, `sphere_hit_shadow`) that are slightly more optimal because they don't need to calculate normals and can work with squared distances.
-   It doesn't check for intersection with the light source itself.

### `initialize_scene()`

-   Defines materials with different colors and types (`LAMBERTIAN`, `DIFFUSE_LIGHT`).
-   Defines the geometry of the scene: 5 walls for the Cornell Box, a light source box, and two spheres.

## 4. Utility Components

### Vector Math (`vec3`, `vec2`)

-   A set of functions for basic vector arithmetic: addition, subtraction, scaling, dot product, cross product, length, and normalization.

### Math Approximations

-   `sqrtf_approx`, `sqrtf_accurate_approx`: Faster square root calculations.
-   `fast_inv_sqrt`: The famous "fast inverse square root" algorithm from Quake III.
-   `fast_inv`: Fast reciprocal (1/x) using the inverse square root trick.

### Random Number Generation

-   `hash2`, `hash3`: Simple hashing functions that produce pseudo-random 2D and 3D vectors based on a float `seed`. This provides deterministic "randomness" needed for reproducible images.
-   `random_cos_weighted_hemisphere_direction`: Generates a random direction on a hemisphere, biased towards the normal, which is important for diffuse reflections.
-   `random_in_unit_disk`: Generates a random point in a disk, used for simulating a camera lens (depth of field).

### Data Structures

-   `ray`: `{origin, direction}`
-   `material`: `{type, color}`
-   `hit_record`: Stores all information about a ray-object intersection.
-   `sphere_primitive`, `box_primitive`: Stores the geometric data for primitives.
-   `camera`: Stores all camera parameters (position, orientation, lens settings). 