#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Path tracing settings
#define SAMPLES_PER_PIXEL 32
#define MAX_BOUNCES 3
#define FOV 60.0

// Scene dimensions
#define NUM_SPHERES 2
#define NUM_PLANES 6

// Fixed-point math settings
#define FRAC_BITS 16
#define ONE (1 << FRAC_BITS)
#define F(x) ((int32_t)((x) * ONE))
#define I(x) ((x) >> FRAC_BITS)

// Simple vector struct
typedef struct {
    int32_t x, y, z;
} Vec3;

// Ray
typedef struct {
    Vec3 orig, dir;
} Ray;

// Material
typedef struct {
    Vec3 color;
    int is_light;
} Material;

// Sphere
typedef struct {
    Vec3 center;
    int32_t radius;
    Material material;
} Sphere;

// Plane
typedef struct {
    Vec3 normal;
    int32_t dist;
    Material material;
} Plane;

// Structure to hold intersection results
typedef struct {
    int32_t t;
    int hit_type; // 0 for sphere, 1 for plane
    int hit_index;
    int hit;      // 1 if an object was hit, 0 otherwise
} Intersection;

#define UNIT_VECTOR_LUT_SIZE 64
static const Vec3 g_unit_vector_lut[UNIT_VECTOR_LUT_SIZE] = {
    {.x = 37813, .y = 37813, .z = 37813}, {.x = -37813, .y = -37813, .z = -37813},
    {.x = 53475, .y = 0, .z = 37813}, {.x = -53475, .y = 0, .z = -37813},
    {.x = 0, .y = 53475, .z = 37813}, {.x = 0, .y = -53475, .z = -37813},
    {.x = 37813, .y = 53475, .z = 0}, {.x = -37813, .y = -53475, .z = 0},
    {.x = 65536, .y = 0, .z = 0}, {.x = -65536, .y = 0, .z = 0},
    {.x = 0, .y = 65536, .z = 0}, {.x = 0, .y = -65536, .z = 0},
    {.x = 0, .y = 0, .z = 65536}, {.x = 0, .y = 0, .z = -65536},
    {.x = 46340, .y = 46340, .z = 0}, {.x = -46340, .y = -46340, .z = 0},
    {.x = 46340, .y = 0, .z = 46340}, {.x = -46340, .y = 0, .z = -46340},
    {.x = 0, .y = 46340, .z = 46340}, {.x = 0, .y = -46340, .z = -46340},
    {.x = 25881, .y = 25881, .z = 57925}, {.x = -25881, .y = -25881, .z = -57925},
    {.x = 25881, .y = 57925, .z = 25881}, {.x = -25881, .y = -57925, .z = -25881},
    {.x = 57925, .y = 25881, .z = 25881}, {.x = -57925, .y = -25881, .z = -25881},
    {.x = 32768, .y = 32768, .z = 46340}, {.x = -32768, .y = -32768, .z = -46340},
    {.x = 32768, .y = 46340, .z = 32768}, {.x = -32768, .y = -46340, .z = -32768},
    {.x = 46340, .y = 32768, .z = 32768}, {.x = -46340, .y = -32768, .z = -32768},
    {.x = 12345, .y = 64543, .z = 10000}, {.x = -12345, .y = -64543, .z = -10000},
    {.x = 64543, .y = 10000, .z = 12345}, {.x = -64543, .y = -10000, .z = -12345},
    {.x = 10000, .y = 12345, .z = 64543}, {.x = -10000, .y = -12345, .z = -64543},
    {.x = 40000, .y = 50000, .z = 1234}, {.x = -40000, .y = -50000, .z = -1234},
    {.x = 50000, .y = 1234, .z = 40000}, {.x = -50000, .y = -1234, .z = -40000},
    {.x = 1234, .y = 40000, .z = 50000}, {.x = -1234, .y = -40000, .z = -50000},
    {.x = 18000, .y = 18000, .z = 60963}, {.x = -18000, .y = -18000, .z = -60963},
    {.x = 18000, .y = 60963, .z = 18000}, {.x = -18000, .y = -60963, .z = -18000},
    {.x = 60963, .y = 18000, .z = 18000}, {.x = -60963, .y = -18000, .z = -18000},
    {.x = 55555, .y = 22222, .z = 25000}, {.x = -55555, .y = -22222, .z = -25000},
    {.x = 22222, .y = 25000, .z = 55555}, {.x = -22222, .y = -25000, .z = -55555},
    {.x = 25000, .y = 55555, .z = 22222}, {.x = -25000, .y = -55555, .z = -22222},
    {.x = 61234, .y = 12345, .z = 20000}, {.x = -61234, .y = -12345, .z = -20000},
    {.x = 12345, .y = 20000, .z = 61234}, {.x = -12345, .y = -20000, .z = -61234},
    {.x = 20000, .y = 61234, .z = 12345}, {.x = -20000, .y = -61234, .z = -12345}
};

static const Sphere g_spheres[NUM_SPHERES] = {
    {.center = {.x = F(-0.6), .y = F(0.0), .z = F(-2.8)}, .radius = F(0.7), .material = {.color = {.x = F(1.0), .y = F(0.7), .z = F(0.2)}, .is_light = 0}}, // Large yellow sphere
    {.center = {.x = F(0.6), .y = F(-0.5), .z = F(-3.2)}, .radius = F(0.5), .material = {.color = {.x = F(0.5), .y = F(0.5), .z = F(0.5)}, .is_light = 0}}  // Small grey sphere
};

static const Plane g_planes[NUM_PLANES] = {
    {.normal = {.x = F(0), .y = F(1), .z = F(0)}, .dist = F(-1), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},   // Floor
    // Emissive light panel – slightly below the ceiling so the ceiling itself can receive light
    {.normal = {.x = F(0), .y = F(-1), .z = F(0)}, .dist = F(-2.99), .material = {.color = {.x = F(32.0), .y = F(32.0), .z = F(32.0)}, .is_light = 1}},
    // Actual ceiling (non-emissive)
    {.normal = {.x = F(0), .y = F(-1), .z = F(0)}, .dist = F(-3), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},
    {.normal = {.x = F(1), .y = F(0), .z = F(0)}, .dist = F(-2), .material = {.color = {.x = F(0.75), .y = F(0.25), .z = F(0.25)}, .is_light = 0}},         // Left wall (red)
    {.normal = {.x = F(-1), .y = F(0), .z = F(0)}, .dist = F(-2), .material = {.color = {.x = F(0.25), .y = F(0.75), .z = F(0.25)}, .is_light = 0}},        // Right wall (green)
    {.normal = {.x = F(0), .y = F(0), .z = F(1)}, .dist = F(-5), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},   // Back wall
};

// Fixed-point multiplication
int32_t mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) >> FRAC_BITS);
}

// Fixed-point division
int32_t div_fp(int32_t a, int32_t b) {
    if (b == 0) return 0;
    //int32_t X0 = 
    return (int32_t)(((int64_t)a << FRAC_BITS) / b);
}


// Fast inverse square root for fixed-point numbers
static inline int32_t inv_sqrt_fp(int32_t x) {
    if (x <= 0) return 0;
    float x_f = (float)x / (float)ONE;
    union { float f; uint32_t i; } u;
    u.f = x_f;
    u.i = 0x5f3759df - (u.i >> 1);
    u.f = u.f * (1.5f - 0.5f * x_f * u.f * u.f);
    return (int32_t)(u.f * (float)ONE);
}

// Fixed-point square root
int32_t sqrt_fp(int32_t n) {
    if (n <= 0) return 0;
    return mul(n, inv_sqrt_fp(n));
}

// Random number generation
static uint32_t rand_state = 12345;
static uint32_t rand_u32() { // maybe make it generate 3 randon numbers each clock cycle?
    // xorshift
    rand_state ^= rand_state << 13;
    rand_state ^= rand_state >> 17;
    rand_state ^= rand_state << 5;
    return rand_state;
}

int32_t rand_fp() {
    return (int32_t)((uint64_t)rand_u32() * ONE >> 32);
}


Vec3 random_unit_vector() {
    uint32_t r_val = rand_u32();
    int lut_idx = r_val % UNIT_VECTOR_LUT_SIZE;
    return g_unit_vector_lut[lut_idx];
}

// Forward declarations
int is_on_light(Vec3 p);
int32_t intersect_sphere(Ray r, Sphere s);
int32_t intersect_plane(Ray r, Plane p);
Intersection intersect_scene(Ray r);
Vec3 trace_path(Ray r);


// Vector operations
Vec3 vec_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vec_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
int32_t vec_dot(Vec3 a, Vec3 b) { return mul(a.x, b.x) + mul(a.y, b.y) + mul(a.z, b.z); }
Vec3 vec_mul(Vec3 a, Vec3 b) { return (Vec3){mul(a.x, b.x), mul(a.y, b.y), mul(a.z, b.z)}; }
Vec3 vec_scale(Vec3 v, int32_t s) { return (Vec3){mul(v.x, s), mul(v.y, s), mul(v.z, s)}; }
int32_t vec_len_sq(Vec3 v) { return vec_dot(v, v); }
Vec3 vec_norm(Vec3 v) {
    int32_t len_sq = vec_len_sq(v);
    int32_t inv_len = inv_sqrt_fp(len_sq);
    return vec_scale(v, inv_len);
}

// Returns an Intersection result.
Intersection intersect_scene(Ray r) {
    Intersection result;
    result.t = F(1e9);
    result.hit_type = -1;
    result.hit_index = -1;
    result.hit = 0;

    // Sphere intersection
    for (size_t i = 0; i < NUM_SPHERES; ++i) {
        int32_t t = intersect_sphere(r, g_spheres[i]);
        if (t < result.t) {
            result.t = t;
            result.hit_type = 0; // 0 for sphere
            result.hit_index = (int)i;
        }
    }
    
    // Plane intersection
    for (size_t i = 0; i < NUM_PLANES; ++i) {
        int32_t t = intersect_plane(r, g_planes[i]);
        if (t < result.t) {
            result.t = t;
            result.hit_type = 1; // 1 for plane
            result.hit_index = (int)i;
        }
    }
    
    if (result.hit_type != -1) {
        result.hit = 1;
    }
    
    return result;
}

// Ray-sphere intersection
int32_t intersect_sphere(Ray r, Sphere s) {
    Vec3 oc = vec_sub(r.orig, s.center);
    int32_t a = vec_dot(r.dir, r.dir);
    int32_t b = 2 * vec_dot(oc, r.dir);
    int32_t c = vec_dot(oc, oc) - mul(s.radius, s.radius);
    int32_t discriminant = mul(b, b) - 4 * mul(a, c);
    if (discriminant < 0) return F(1e9);
    
    int32_t sqrt_d = sqrt_fp(discriminant);
    // Pre-calculate inverse of denominator to replace two divisions with one.
    int32_t inv_2a = div_fp(ONE, 2 * a);

    int32_t t = mul(-b - sqrt_d, inv_2a);
    if (t > F(0.001)) return t;
    
    int32_t t2 = mul(-b + sqrt_d, inv_2a);
    if (t2 > F(0.001)) return t2;
    
    return F(1e9);
}

// Ray-plane intersection
int32_t intersect_plane(Ray r, Plane p) {
    int32_t denom = vec_dot(p.normal, r.dir);
    if (denom > -F(0.001) && denom < F(0.001)) return F(1e9); // Parallel
    int32_t t = div_fp(vec_dot(p.normal, vec_sub(vec_scale(p.normal, p.dist), r.orig)), denom);
    if (t <= F(0.001)) return F(1e9);

    if (p.material.is_light) {
        Vec3 hit_pt = vec_add(r.orig, vec_scale(r.dir, t));
        if (!is_on_light(hit_pt)) {
            return F(1e9);
        }
    }

    return t;
}

Vec3 trace_path(Ray r) {

    Vec3 path_color = {F(0), F(0), F(0)};
    Vec3 path_attenuation = {ONE, ONE, ONE};

    for (int b = 0; b < MAX_BOUNCES; ++b) {
        Intersection inter = intersect_scene(r);

        if (!inter.hit) {
            path_attenuation = (Vec3){F(0), F(0), F(0)};
        }

        int32_t t = inter.t;
        int hit_object_type = inter.hit_type;
        int hit_object_index = inter.hit_index;

        Vec3 hit_point = vec_add(r.orig, vec_scale(r.dir, t));
        Vec3 hit_normal;
        Material mat;

        if (hit_object_type == 0) { // Sphere
            mat = g_spheres[hit_object_index].material;
            hit_normal = vec_norm(vec_sub(hit_point, g_spheres[hit_object_index].center));
        } else { // Plane
            mat = g_planes[hit_object_index].material;
            hit_normal = g_planes[hit_object_index].normal;
        }

        Material surface_mat = mat;
        if (mat.is_light) { // If we hit the ceiling plane
            if (is_on_light(hit_point)) {
                // Only add emission for camera rays to avoid double counting with NEE
                if (b == 0) {
                    path_color = vec_add(path_color, mat.color);
                }
                path_attenuation = (Vec3){F(0), F(0), F(0)};
            } else {
                // Hit ceiling, but outside the light. Treat as grey.
                surface_mat.color = (Vec3){F(0.2), F(0.2), F(0.2)};
            }
        }

        // Check if it is in a shadow
        uint32_t r_val = rand_u32();
        int32_t rand1 = r_val & 0xFFFF;
        int32_t rand2 = r_val >> 16;
        Vec3 light_point = {F(-1.0) + mul(F(2.0), rand1), F(2.99), F(-3.2) + mul(F(0.4), rand2)};
        Vec3 light_vec = vec_sub(light_point, hit_point);
        int32_t dist_sq = vec_len_sq(light_vec);
        Vec3 light_dir = vec_norm(light_vec);

        Ray shadow_ray = {vec_add(hit_point, vec_scale(hit_normal, F(0.01))), light_dir};
        int occluded = 0;
        for (size_t i = 0; i < NUM_SPHERES; ++i) {
            int32_t shadow_t = intersect_sphere(shadow_ray, g_spheres[i]);
            if (shadow_t < F(1e8) && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; break; }
        }
        for (size_t i = 0; i < NUM_PLANES; ++i) {
            if (g_planes[i].material.is_light) continue; // Don't treat the emissive plane as occluder
            int32_t shadow_t = intersect_plane(shadow_ray, g_planes[i]);
            if (shadow_t < F(1e8) && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; break; }
        }

        if (!occluded) {
            // if it is NOT in a shadow, calculate the direct light contribution
            int32_t cos_theta = vec_dot(hit_normal, light_dir);
            Vec3 light_normal = {0, -ONE, 0}; // Light surface normal points down into the box
            int32_t cos_alpha = vec_dot(light_normal, vec_scale(light_dir, -ONE));

            if (cos_theta > 0 && cos_alpha > 0) { // if the light and surface are facing each other
                Material light_mat = g_planes[1].material;
                int32_t light_area = F(2.0 * 0.4);
                int32_t geom_term_num = mul(cos_theta, cos_alpha);
                int32_t geom_term = div_fp(geom_term_num, dist_sq);

                Vec3 direct_light = vec_mul(path_attenuation, surface_mat.color);
                direct_light = vec_mul(direct_light, light_mat.color);
                direct_light = vec_scale(direct_light, geom_term);
                direct_light = vec_scale(direct_light, light_area);
                // Divide by PI for diffuse BRDF
                direct_light = vec_scale(direct_light, F(0.3183)); 
                path_color = vec_add(path_color, direct_light);
            }
        }

        // Attenuate path for next bounce (indirect light)
        path_attenuation = vec_mul(path_attenuation, surface_mat.color);
        
        // New random direction for bounced ray
        Vec3 random_dir = random_unit_vector();
        Vec3 bounce_dir = vec_add(hit_normal, random_dir);

        r.orig = vec_add(hit_point, vec_scale(hit_normal, F(0.01)));
        r.dir = vec_norm(bounce_dir);
    }

    return path_color;
}

// Check if a point is on the rectangular light source on the ceiling
int is_on_light(Vec3 p) {
    return (p.x >= F(-1.0) && p.x <= F(1.0) && p.z >= F(-3.2) && p.z <= F(-2.8));
}


// EVERYTHING BELOW HERE WILL BE IMPLEMENTED IN THE HOST MACHINE

// PPM image output
void write_ppm(const char *filename, int width, int height, const uint8_t *data) {
    FILE *f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(data, 1, width * height * 3, f);
    fclose(f);
}

int main() {
    int width = 512;
    int height = 512;
    uint8_t *image = (uint8_t *)malloc(width * height * 3);

    // Scene is now defined globally and is constant.

    // Camera
    Ray cam = {{F(0), F(0.8), F(2)}, {F(0), F(0), F(-1)}};
    double fov_rad = FOV * M_PI / 180.0;
    double fov_scale = tan(fov_rad / 2.0);
    double aspect_ratio = (double)width / height;

    for (int y = 0; y < height; ++y) {
        fprintf(stderr, "\rRendering: %3d%%", (y * 100) / (height - 1));
        fflush(stderr);
        for (int x = 0; x < width; ++x) {
            Vec3 pixel_color = {F(0), F(0), F(0)};

            for (int s = 0; s < SAMPLES_PER_PIXEL; ++s) {
                // Screen coordinates with antialiasing
                double jitter_x = (double)rand_fp() / ONE;
                double jitter_y = (double)rand_fp() / ONE;

                double sx_ndc = (2.0 * (x + jitter_x) / width) - 1.0;
                double sy_ndc = 1.0 - (2.0 * (y + jitter_y) / height);
                
                Ray r = cam;
                r.dir.x = F(sx_ndc * aspect_ratio * fov_scale);
                r.dir.y = F(sy_ndc * fov_scale);
                r.dir = vec_norm(r.dir);

                pixel_color = vec_add(pixel_color, trace_path(r));
            }
            
            Vec3 color = vec_scale(pixel_color, div_fp(ONE, F(SAMPLES_PER_PIXEL)));

            int i = (y * width + x) * 3;
            int r_val = ((int64_t)color.x * 255) >> FRAC_BITS;
            int g_val = ((int64_t)color.y * 255) >> FRAC_BITS;
            int b_val = ((int64_t)color.z * 255) >> FRAC_BITS;
            
            image[i]   = r_val > 255 ? 255 : (r_val < 0 ? 0 : r_val);
            image[i+1] = g_val > 255 ? 255 : (g_val < 0 ? 0 : g_val);
            image[i+2] = b_val > 255 ? 255 : (b_val < 0 ? 0 : b_val);
        }
    }
    fprintf(stderr, "\nDone.\n");

    write_ppm("output.ppm", width, height, image);
    free(image);

    return 0;
} 