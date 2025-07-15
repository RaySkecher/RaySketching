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
#define FOV 60.0 // Field of View in degrees

// Fixed-point math settings
#define FRAC_BITS 16
#define ONE (1 << FRAC_BITS)
#define F(x) ((int32_t)((x) * ONE))
#define I(x) ((x) >> FRAC_BITS)

// Lookup table for random angles
#define ANGLE_LUT_SIZE 256
static int32_t g_cos_lut[ANGLE_LUT_SIZE];
static int32_t g_sin_lut[ANGLE_LUT_SIZE];

void init_angle_lut() {
    for (int i = 0; i < ANGLE_LUT_SIZE; ++i) {
        double angle = 2.0 * M_PI * i / ANGLE_LUT_SIZE;
        g_cos_lut[i] = F(cos(angle));
        g_sin_lut[i] = F(sin(angle));
    }
}

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

// Fixed-point multiplication
int32_t mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) >> FRAC_BITS);
}

// Fixed-point division
int32_t div_fp(int32_t a, int32_t b) {
    if (b == 0) return 0;
    return (int32_t)(((int64_t)a << FRAC_BITS) / b);
}

// Random number generation
static uint32_t rand_state = 12345;
static uint32_t rand_u32() {
    // xorshift
    rand_state ^= rand_state << 13;
    rand_state ^= rand_state >> 17;
    rand_state ^= rand_state << 5;
    return rand_state;
}
int32_t rand_fp() {
    return (int32_t)((uint64_t)rand_u32() * ONE >> 32);
}

// Forward declarations
Vec3 random_unit_vector();
int is_on_light(const Vec3* p);
int32_t intersect_sphere(const Ray *r, const Sphere *s);
int32_t intersect_plane(const Ray *r, const Plane *p);

// Fixed-point square root
int32_t sqrt_fp(int32_t n) {
    if (n <= 0) {
        return 0;
    }

    int32_t root = n;
    if (root < ONE) {
        root = ONE;
    }

    for (int i = 0; i < 20; ++i) {
        int32_t last_root = root;
        if (last_root == 0) { // Should not happen for n > 0
            break;
        }
        root = (root + div_fp(n, root)) >> 1;
        if (root == last_root) {
            break;
        }
    }

    return root;
}

// Vector operations
Vec3 vec_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vec_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
int32_t vec_dot(Vec3 a, Vec3 b) { return mul(a.x, b.x) + mul(a.y, b.y) + mul(a.z, b.z); }
Vec3 vec_mul(Vec3 a, Vec3 b) { return (Vec3){mul(a.x, b.x), mul(a.y, b.y), mul(a.z, b.z)}; }
Vec3 vec_scale(Vec3 v, int32_t s) { return (Vec3){mul(v.x, s), mul(v.y, s), mul(v.z, s)}; }
int32_t vec_len_sq(Vec3 v) { return vec_dot(v, v); }
Vec3 vec_norm(Vec3 v) {
    int32_t len_sq = vec_len_sq(v);
    if (len_sq == 0) return v;
    int32_t len = sqrt_fp(len_sq);
    if (len == 0) return v;
    return vec_scale(v, div_fp(ONE, len));
}

// Returns 1 if an object was hit, 0 otherwise.
int intersect_scene(const Ray *r,
                    const Sphere* spheres, size_t num_spheres,
                    const Plane* planes, size_t num_planes,
                    int32_t* t_out, int* hit_type_out, int* hit_index_out) {
    
    int32_t min_t = F(1e9);
    *hit_type_out = -1;
    *hit_index_out = -1;

    // Sphere intersection
    for (size_t i = 0; i < num_spheres; ++i) {
        int32_t t = intersect_sphere(r, &spheres[i]);
        if (t < min_t) {
            min_t = t;
            *hit_type_out = 0; // 0 for sphere
            *hit_index_out = (int)i;
        }
    }
    
    // Plane intersection
    for (size_t i = 0; i < num_planes; ++i) {
        int32_t t = intersect_plane(r, &planes[i]);
        if (t < min_t) {
            min_t = t;
            *hit_type_out = 1; // 1 for plane
            *hit_index_out = (int)i;
        }
    }
    
    if (*hit_type_out != -1) {
        *t_out = min_t;
        return 1;
    }
    
    return 0;
}

// Ray-sphere intersection
int32_t intersect_sphere(const Ray *r, const Sphere *s) {
    Vec3 oc = vec_sub(r->orig, s->center);
    int32_t a = vec_dot(r->dir, r->dir);
    int32_t b = 2 * vec_dot(oc, r->dir);
    int32_t c = vec_dot(oc, oc) - mul(s->radius, s->radius);
    int32_t discriminant = mul(b, b) - 4 * mul(a, c);
    if (discriminant < 0) return F(1e9);
    
    int32_t t = div_fp(-b - sqrt_fp(discriminant), 2 * a);
    if (t > F(0.001)) return t;
    
    t = div_fp(-b + sqrt_fp(discriminant), 2 * a);
    if (t > F(0.001)) return t;
    
    return F(1e9);
}

// Ray-plane intersection
int32_t intersect_plane(const Ray *r, const Plane *p) {
    int32_t denom = vec_dot(p->normal, r->dir);
    if (denom > -F(0.001) && denom < F(0.001)) return F(1e9); // Parallel
    int32_t t = div_fp(vec_dot(p->normal, vec_sub(vec_scale(p->normal, p->dist), r->orig)), denom);
    if (t <= F(0.001)) return F(1e9);

    if (p->material.is_light) {
        Vec3 hit_pt = vec_add(r->orig, vec_scale(r->dir, t));
        if (!is_on_light(&hit_pt)) {
            return F(1e9);
        }
    }

    return t;
}

Vec3 trace_path(Ray r,
                const Sphere* spheres, size_t num_spheres,
                const Plane* planes, size_t num_planes) {

    Vec3 path_color = {F(0), F(0), F(0)};
    Vec3 path_attenuation = {ONE, ONE, ONE};

    for (int b = 0; b < MAX_BOUNCES; ++b) {
        int32_t t;
        int hit_object_type;
        int hit_object_index;

        if (!intersect_scene(&r, spheres, num_spheres, planes, num_planes, &t, &hit_object_type, &hit_object_index)) {
            break; // Ray escaped
        }

        Vec3 hit_point = vec_add(r.orig, vec_scale(r.dir, t));
        Vec3 hit_normal;
        Material mat;

        if (hit_object_type == 0) { // Sphere
            mat = spheres[hit_object_index].material;
            hit_normal = vec_norm(vec_sub(hit_point, spheres[hit_object_index].center));
        } else { // Plane
            mat = planes[hit_object_index].material;
            hit_normal = planes[hit_object_index].normal;
        }

        Material surface_mat = mat;
        if (mat.is_light) { // If we hit the ceiling plane
            if (is_on_light(&hit_point)) {
                // Only add emission for camera rays to avoid double counting with NEE
                if (b == 0) {
                    path_color = vec_add(path_color, mat.color);
                }
                break; // Path ends if it hits a light
            } else {
                // Hit ceiling, but outside the light. Treat as grey.
                surface_mat.color = (Vec3){F(0.2), F(0.2), F(0.2)};
            }
        }

        // Next Event Estimation (direct light sampling)
        {
            Vec3 light_point = {F(-1.0) + mul(F(2.0), rand_fp()), F(2.99), F(-3.2) + mul(F(0.4), rand_fp())};
            Vec3 light_vec = vec_sub(light_point, hit_point);
            int32_t dist_sq = vec_len_sq(light_vec);
            Vec3 light_dir = vec_norm(light_vec);

            Ray shadow_ray = {vec_add(hit_point, vec_scale(hit_normal, F(0.001))), light_dir};
            int occluded = 0;
            for (size_t i = 0; i < num_spheres; ++i) {
                int32_t shadow_t = intersect_sphere(&shadow_ray, &spheres[i]);
                if (shadow_t < F(1e8) && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; break; }
            }
            if (!occluded) {
                for (size_t i = 0; i < num_planes; ++i) {
                    if (planes[i].material.is_light) continue; // Don't treat the emissive plane as occluder
                    int32_t shadow_t = intersect_plane(&shadow_ray, &planes[i]);
                    if (shadow_t < F(1e8) && mul(shadow_t, shadow_t) < dist_sq) { occluded = 1; break; }
                }
            }

            if (!occluded) {
                int32_t cos_theta = vec_dot(hit_normal, light_dir);
                Vec3 light_normal = {0, -ONE, 0}; // Light surface normal points down into the box
                int32_t cos_alpha = vec_dot(light_normal, vec_scale(light_dir, -ONE));

                if (cos_theta > 0 && cos_alpha > 0) {
                    Material light_mat = planes[1].material;
                    int32_t light_area = F(2.0 * 0.4);
                    int32_t geom_term_num = mul(cos_theta, cos_alpha);
                    int32_t geom_term = div_fp(geom_term_num, dist_sq);

                    if (geom_term > 0) {
                        Vec3 direct_light = vec_mul(path_attenuation, surface_mat.color);
                        direct_light = vec_mul(direct_light, light_mat.color);
                        direct_light = vec_scale(direct_light, geom_term);
                        direct_light = vec_scale(direct_light, light_area);
                        // Divide by PI for diffuse BRDF
                        direct_light = vec_scale(direct_light, F(0.3183)); 
                        path_color = vec_add(path_color, direct_light);
                    }
                }
            }
        }

        // Attenuate path for next bounce (indirect light)
        path_attenuation = vec_mul(path_attenuation, surface_mat.color);
        
        // New random direction for bounced ray
        Vec3 random_dir = random_unit_vector();
        Vec3 bounce_dir = vec_add(hit_normal, random_dir);
        if (vec_len_sq(bounce_dir) == 0) {
            bounce_dir = hit_normal;
        }

        r.orig = vec_add(hit_point, vec_scale(hit_normal, F(0.001)));
        r.dir = vec_norm(bounce_dir);
    }

    return path_color;
}

// Check if a point is on the rectangular light source on the ceiling
int is_on_light(const Vec3* p) {
    return (p->x >= F(-1.0) && p->x <= F(1.0) && p->z >= F(-3.2) && p->z <= F(-2.8));
}

Vec3 random_unit_vector() {
    // Generate a random angle phi
    uint32_t r_val = rand_u32();
    int lut_idx = r_val % ANGLE_LUT_SIZE;
    int32_t cos_phi = g_cos_lut[lut_idx];
    int32_t sin_phi = g_sin_lut[lut_idx];

    // Generate a random z-coordinate (cos_theta)
    int32_t cos_theta = (2 * rand_fp()) - ONE; // Uniform in [-1, 1]
    int32_t sin_theta = sqrt_fp(ONE - mul(cos_theta, cos_theta));

    // Construct the unit vector
    Vec3 p;
    p.x = mul(sin_theta, cos_phi);
    p.y = mul(sin_theta, sin_phi);
    p.z = cos_theta;
    return p;
}

// PPM image output
void write_ppm(const char *filename, int width, int height, const uint8_t *data) {
    FILE *f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(data, 1, width * height * 3, f);
    fclose(f);
}

int main() {
    init_angle_lut();
    int width = 512;
    int height = 512;
    uint8_t *image = (uint8_t *)malloc(width * height * 3);

    // Scene
    Sphere spheres[] = {
        {{F(-0.6), F(0.0), F(-2.8)}, F(0.7), {{F(1.0), F(0.7), F(0.2)}, 0}}, // Large yellow sphere
        {{F(0.6), F(-0.5), F(-3.2)}, F(0.5), {{F(0.5), F(0.5), F(0.5)}, 0}}  // Small grey sphere
    };
    Plane planes[] = {
        {{F(0), F(1), F(0)}, F(-1), {{F(0.75), F(0.75), F(0.75)}, 0}},   // Floor
        // Emissive light panel – slightly below the ceiling so the ceiling itself can receive light
        {{F(0), F(-1), F(0)}, F(-2.99), {{F(32.0), F(32.0), F(32.0)}, 1}},
        // Actual ceiling (non-emissive)
        {{F(0), F(-1), F(0)}, F(-3),    {{F(0.75), F(0.75), F(0.75)}, 0}},
        {{F(1), F(0), F(0)}, F(-2), {{F(0.75), F(0.25), F(0.25)}, 0}},         // Left wall (red)
        {{F(-1), F(0), F(0)}, F(-2), {{F(0.25), F(0.75), F(0.25)}, 0}},        // Right wall (green)
        {{F(0), F(0), F(1)}, F(-5), {{F(0.75), F(0.75), F(0.75)}, 0}},   // Back wall
    };
    const size_t num_spheres = sizeof(spheres) / sizeof(Sphere);
    const size_t num_planes = sizeof(planes) / sizeof(Plane);

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

                pixel_color = vec_add(pixel_color, trace_path(r, spheres, num_spheres, planes, num_planes));
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