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
#define LIGHT_INTENSITY 2.0

// Scene dimensions
#define NUM_SPHERES 2
#define NUM_PLANES 6

// Fixed-point math settings — 16-bit total (4 integer + 12 fractional)
typedef int16_t fp_t;
#define FRAC_BITS 12
#define ONE (1 << FRAC_BITS)
// Convert a floating-point value to 4.12 fixed-point (rounded)
#define F(x) ((fp_t)((x) * ONE + ((x) >= 0 ? 0.5f : -0.5f)))
#define I(x) ((x) >> FRAC_BITS)

// Handy constants
#define FP_EPS ((fp_t)1)             // ≈ 0.00024 in real units
#define FP_INF 0x7FFFFFFF            // Large “infinite” distance sentinel

// Simple vector struct
typedef struct {
    fp_t x, y, z;
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
    fp_t radius;
    Material material;
} Sphere;

// Plane
typedef struct {
    Vec3 normal;
    fp_t dist;
    Material material;
} Plane;

// Structure to hold intersection results
typedef struct {
    int32_t t;            // keep 32-bit for extra head-room during comparisons
    int hit_type; // 0 for sphere, 1 for plane
    int hit_index;
    int hit;      // 1 if an object was hit, 0 otherwise
} Intersection;

#define UNIT_VECTOR_LUT_SIZE 128
// 64 pre-computed unit vectors stored at 8.8 precision; they will be left-shifted
// at runtime to 4.12 to keep the table small and readable.
static const Vec3 g_unit_vector_lut[UNIT_VECTOR_LUT_SIZE] = {
    {.x = 148,  .y = 148,  .z = 148 }, {.x = -148, .y = -148, .z = -148},
    {.x = 210,  .y = 0,    .z = 148 }, {.x = -210, .y = 0,    .z = -148},
    {.x = 0,    .y = 210,  .z = 148 }, {.x = 0,    .y = -210, .z = -148},
    {.x = 148,  .y = 210,  .z = 0   }, {.x = -148, .y = -210, .z = 0   },
    {.x = 256,  .y = 0,    .z = 0   }, {.x = -256, .y = 0,    .z = 0   },
    {.x = 0,    .y = 256,  .z = 0   }, {.x = 0,    .y = -256, .z = 0   },
    {.x = 0,    .y = 0,    .z = 256 }, {.x = 0,    .y = 0,    .z = -256},
    {.x = 182,  .y = 182,  .z = 0   }, {.x = -182, .y = -182, .z = 0   },
    {.x = 182,  .y = 0,    .z = 182 }, {.x = -182, .y = 0,    .z = -182},
    {.x = 0,    .y = 182,  .z = 182 }, {.x = 0,    .y = -182, .z = -182},
    {.x = 102,  .y = 102,  .z = 226 }, {.x = -102, .y = -102, .z = -226},
    {.x = 102,  .y = 226,  .z = 102 }, {.x = -102, .y = -226, .z = -102},
    {.x = 226,  .y = 102,  .z = 102 }, {.x = -226, .y = -102, .z = -102},
    {.x = 128,  .y = 128,  .z = 182 }, {.x = -128, .y = -128, .z = -182},
    {.x = 128,  .y = 182,  .z = 128 }, {.x = -128, .y = -182, .z = -128},
    {.x = 182,  .y = 128,  .z = 128 }, {.x = -182, .y = -128, .z = -128},
    {.x =  48,  .y = 252,  .z =  39 }, {.x =  -48, .y = -252, .z =  -39},
    {.x = 252,  .y =  39,  .z =  48 }, {.x = -252, .y =  -39, .z =  -48},
    {.x =  39,  .y =  48,  .z = 252 }, {.x =  -39, .y =  -48, .z = -252},
    {.x = 156,  .y = 195,  .z =   5 }, {.x = -156, .y = -195, .z =   -5},
    {.x = 195,  .y =   5,  .z = 156 }, {.x = -195, .y =   -5, .z = -156},
    {.x =   5,  .y = 156,  .z = 195 }, {.x =   -5, .y = -156, .z = -195},
    {.x =  70,  .y =  70,  .z = 238 }, {.x =  -70, .y =  -70, .z = -238},
    {.x =  70,  .y = 238,  .z =  70 }, {.x =  -70, .y = -238, .z =  -70},
    {.x = 238,  .y =  70,  .z =  70 }, {.x = -238, .y =  -70, .z =  -70},
    {.x = 217,  .y =  86,  .z =  98 }, {.x = -217, .y =  -86, .z =  -98},
    {.x =  86,  .y =  98,  .z = 217 }, {.x =  -86, .y =  -98, .z = -217},
    {.x =  98,  .y = 217,  .z =  86 }, {.x =  -98, .y = -217, .z =  -86},
    {.x = 239,  .y =  48,  .z =  78 }, {.x = -239, .y =  -48, .z =  -78},
    {.x =  48,  .y =  78,  .z = 239 }, {.x =  -48, .y =  -78, .z = -239},
    {.x =  78,  .y = 239,  .z =  48 }, {.x =  -78, .y = -239, .z =  -48},
    {.x = 180, .y =  34, .z = 193}, {.x = -180, .y =  34, .z = -193},
    {.x =  34, .y = 193, .z = 180}, {.x =  -34, .y = -193, .z = -180},
    {.x = 193, .y = 180, .z =  34}, {.x = -193, .y = -180, .z =  -34},
    {.x =  45, .y = 210, .z = 160}, {.x =  -45, .y = -210, .z = -160},
    {.x = 210, .y = 160, .z =  45}, {.x = -210, .y = -160, .z =  -45},
    {.x = 160, .y =  45, .z = 210}, {.x = -160, .y =  -45, .z = -210},
    {.x = 125, .y = 230, .z =  70}, {.x = -125, .y = -230, .z =  -70},
    {.x = 230, .y =  70, .z = 125}, {.x = -230, .y =  -70, .z = -125},
    {.x =  70, .y = 125, .z = 230}, {.x =  -70, .y = -125, .z = -230},
    {.x = 190, .y = 200, .z =  20}, {.x = -190, .y = -200, .z =  -20},
    {.x = 200, .y =  20, .z = 190}, {.x = -200, .y =  -20, .z = -190},
    {.x =  20, .y = 190, .z = 200}, {.x =  -20, .y = -190, .z = -200},
    {.x =  95, .y = 215, .z = 140}, {.x =  -95, .y = -215, .z = -140},
    {.x = 215, .y = 140, .z =  95}, {.x = -215, .y = -140, .z =  -95},
    {.x = 140, .y =  95, .z = 215}, {.x = -140, .y =  -95, .z = -215},
    {.x =  60, .y = 180, .z = 230}, {.x =  -60, .y = -180, .z = -230},
    {.x = 180, .y = 230, .z =  60}, {.x = -180, .y = -230, .z =  -60},
    {.x = 230, .y =  60, .z = 180}, {.x = -230, .y =  -60, .z = -180},
    {.x = 130, .y = 200, .z = 150}, {.x = -130, .y = -200, .z = -150},
    {.x = 200, .y = 150, .z = 130}, {.x = -200, .y = -150, .z = -130},
    {.x = 150, .y = 130, .z = 200}, {.x = -150, .y = -130, .z = -200},
    {.x =  85, .y = 145, .z = 240}, {.x =  -85, .y = -145, .z = -240},
    {.x = 145, .y = 240, .z =  85}, {.x = -145, .y = -240, .z =  -85},
    {.x = 240, .y =  85, .z = 145}, {.x = -240, .y =  -85, .z = -145},
    {.x = 170, .y = 100, .z = 220}, {.x = -170, .y = -100, .z = -220},
    {.x = 100, .y = 220, .z = 170}, {.x = -100, .y = -220, .z = -170},
    {.x = 220, .y = 170, .z = 100}, {.x = -220, .y = -170, .z = -100},
    {.x =  66, .y = 210, .z = 195}, {.x =  -66, .y = -210, .z = -195},
    {.x = 210, .y = 195, .z =  66}, {.x = -210, .y = -195, .z =  -66},
    {.x = 195, .y =  66, .z = 210}, {.x = -195, .y =  -66, .z = -210},
    {.x = 110, .y = 180, .z = 210}, {.x = -110, .y = -180, .z = -210}
};

static const Sphere g_spheres[NUM_SPHERES] = {
    {.center = {.x = F(-0.6), .y = F(0.0), .z = F(-2.8)}, .radius = F(0.7), .material = {.color = {.x = F(1.0), .y = F(0.7), .z = F(0.2)}, .is_light = 0}}, // Large yellow sphere
    {.center = {.x = F(0.6), .y = F(-0.5), .z = F(-3.2)}, .radius = F(0.5), .material = {.color = {.x = F(0.5), .y = F(0.5), .z = F(0.5)}, .is_light = 0}}  // Small grey sphere
};

static const Plane g_planes[NUM_PLANES] = {
    {.normal = {.x = F(0), .y = F(1), .z = F(0)}, .dist = F(-1), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},   // Floor
    // Emissive light panel – slightly below the ceiling so the ceiling itself can receive light
    {.normal = {.x = F(0), .y = F(-1), .z = F(0)}, .dist = F(-2.99), .material = {.color = {.x = F(LIGHT_INTENSITY), .y = F(LIGHT_INTENSITY), .z = F(LIGHT_INTENSITY)}, .is_light = 1}},
    // Actual ceiling (non-emissive)
    {.normal = {.x = F(0), .y = F(-1), .z = F(0)}, .dist = F(-3), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},
    {.normal = {.x = F(1), .y = F(0), .z = F(0)}, .dist = F(-2), .material = {.color = {.x = F(0.75), .y = F(0.25), .z = F(0.25)}, .is_light = 0}},         // Left wall (red)
    {.normal = {.x = F(-1), .y = F(0), .z = F(0)}, .dist = F(-2), .material = {.color = {.x = F(0.25), .y = F(0.75), .z = F(0.25)}, .is_light = 0}},        // Right wall (green)
    {.normal = {.x = F(0), .y = F(0), .z = F(1)}, .dist = F(-5), .material = {.color = {.x = F(0.75), .y = F(0.75), .z = F(0.75)}, .is_light = 0}},   // Back wall
};

// Fixed-point multiplication
// 8.8 fixed-point multiply (fp_t × fp_t → 32-bit)
static inline int32_t mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * (int64_t)b) >> FRAC_BITS);
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
    Vec3 base = g_unit_vector_lut[lut_idx];
    // Convert from 8.8 → 4.12 by left-shifting 4 bits.
    return (Vec3){ (fp_t)(base.x << (FRAC_BITS - 8)),
                   (fp_t)(base.y << (FRAC_BITS - 8)),
                   (fp_t)(base.z << (FRAC_BITS - 8)) };
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
Vec3 vec_mul(Vec3 a, Vec3 b) { return (Vec3){ (fp_t)mul(a.x, b.x), (fp_t)mul(a.y, b.y), (fp_t)mul(a.z, b.z)}; }
Vec3 vec_scale(Vec3 v, int32_t s) {
    fp_t s_fp = (fp_t)s; // clamp / cast the scalar into fp_t range
    return (Vec3){ (fp_t)mul(v.x, s_fp), (fp_t)mul(v.y, s_fp), (fp_t)mul(v.z, s_fp)};
}
int32_t vec_len_sq(Vec3 v) { return vec_dot(v, v); }
Vec3 vec_norm(Vec3 v) {
    int32_t len_sq = vec_len_sq(v);
    int32_t inv_len = inv_sqrt_fp(len_sq);
    return vec_scale(v, inv_len);
}

// Returns an Intersection result.
Intersection intersect_scene(Ray r) {
    Intersection result;
    result.t = FP_INF;
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
    if (discriminant < 0) return FP_INF;
    
    int32_t sqrt_d = sqrt_fp(discriminant);
    // Pre-compute 1/(2a) once.
    int32_t inv_2a = div_fp(ONE, 2 * a);

    int32_t t  = mul(-b - sqrt_d, inv_2a);
    if (t > FP_EPS) return t;

    int32_t t2 = mul(-b + sqrt_d, inv_2a);
    if (t2 > FP_EPS) return t2;
    
    return FP_INF;
}

// Ray-plane intersection
int32_t intersect_plane(Ray r, Plane p) {
    int32_t denom = vec_dot(p.normal, r.dir);
    if (denom > -FP_EPS && denom < FP_EPS) return FP_INF; // Parallel
    int32_t t = div_fp(vec_dot(p.normal, vec_sub(vec_scale(p.normal, p.dist), r.orig)), denom);
    if (t <= FP_EPS) return FP_INF;

    if (p.material.is_light) {
        Vec3 hit_pt = vec_add(r.orig, vec_scale(r.dir, t));
        if (!is_on_light(hit_pt)) {
            return FP_INF;
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
        int32_t rand1 = rand_fp(); // 0…ONE
        int32_t rand2 = rand_fp();
        Vec3 light_point = {F(-1.0) + mul(F(2.0), rand1), F(2.99), F(-3.2) + mul(F(0.4), rand2)};
        Vec3 light_vec = vec_sub(light_point, hit_point);
        int32_t dist_sq = vec_len_sq(light_vec);
        Vec3 light_dir = vec_norm(light_vec);

        Ray shadow_ray = {vec_add(hit_point, vec_scale(hit_normal, F(0.01))), light_dir};
        int occluded = 0;
        for (size_t i = 0; i < NUM_SPHERES; ++i) {
            int32_t shadow_t = intersect_sphere(shadow_ray, g_spheres[i]);
            if (shadow_t < FP_INF && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; }
        }
        for (size_t i = 0; i < NUM_PLANES; ++i) {
            if (g_planes[i].material.is_light) continue; // Don't treat the emissive plane as occluder
            int32_t shadow_t = intersect_plane(shadow_ray, g_planes[i]);
            if (shadow_t < FP_INF && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; }
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

#define BRIGHTNESS_SHIFT 4  // final display multiplier (2^n) to compensate for lower 4.12 range

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
            int32_t acc_r = 0, acc_g = 0, acc_b = 0; // 32-bit accumulation to prevent overflow

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

                Vec3 sample_col = trace_path(r);
                acc_r += sample_col.x;
                acc_g += sample_col.y;
                acc_b += sample_col.z;
            }
            
            // Average and convert to 4.12 fixed-point again
            fp_t color_x = (fp_t)(acc_r / SAMPLES_PER_PIXEL);
            fp_t color_y = (fp_t)(acc_g / SAMPLES_PER_PIXEL);
            fp_t color_z = (fp_t)(acc_b / SAMPLES_PER_PIXEL);

            int32_t disp_r = color_x << BRIGHTNESS_SHIFT;
            int32_t disp_g = color_y << BRIGHTNESS_SHIFT;
            int32_t disp_b = color_z << BRIGHTNESS_SHIFT;
 
            int i = (y * width + x) * 3;
            int r_val = ((int64_t)disp_r * 255) >> FRAC_BITS;
            int g_val = ((int64_t)disp_g * 255) >> FRAC_BITS;
            int b_val = ((int64_t)disp_b * 255) >> FRAC_BITS;
            
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