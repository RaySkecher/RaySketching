#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

// Path tracing settings
#define SAMPLES_PER_PIXEL 32
#define MAX_BOUNCES 5

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
int32_t rand_fp() {
    // xorshift
    rand_state ^= rand_state << 13;
    rand_state ^= rand_state >> 17;
    rand_state ^= rand_state << 5;
    return (int32_t)((uint64_t)rand_state * ONE >> 32);
}


// Fixed-point square root
int32_t sqrt_fp(int32_t n) {
    if (n <= 0) {
        return 0;
    }

    int32_t root = n;
    // For values < 1.0, n is a bad initial guess because sqrt(n) > n.
    // Using 1.0 as a starting point is more stable.
    if (root < ONE) {
        root = ONE;
    }

    // Use a bounded number of iterations to prevent infinite loops.
    // Newton's method converges quadratically, so this is plenty.
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

Vec3 random_unit_vector() {
    while (1) {
        Vec3 p = {2 * rand_fp() - ONE, 2 * rand_fp() - ONE, 2 * rand_fp() - ONE};
        int32_t len_sq = vec_len_sq(p);
        if (len_sq > F(0.001) && len_sq < ONE) {
            return vec_norm(p);
        }
    }
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
    return t > F(0.001) ? t : F(1e9);
}

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

    // Scene
    Sphere spheres[] = {
        {{F(-0.5), F(-0.5), F(-3)}, F(0.5), {{F(0.9), F(0.9), F(0.9)}, 0}}, // White sphere
        {{F(0.5), F(0), F(-2.5)}, F(0.7), {{F(0.5), F(0.5), F(1)}, 0}}      // Light blue sphere
    };
    Plane planes[] = {
        {{F(0), F(1), F(0)}, F(-1), {{F(0.8), F(0.8), F(0.8)}, 0}},   // Floor
        {{F(0), F(-1), F(0)}, F(-3), {{F(1.0), F(1.0), F(0.9)}, 1}},  // Ceiling Light
        {{F(1), F(0), F(0)}, F(-2), {{F(1), F(0), F(0)}, 0}},         // Left wall (red)
        {{F(-1), F(0), F(0)}, F(-2), {{F(0), F(1), F(0)}, 0}},        // Right wall (green)
        {{F(0), F(0), F(1)}, F(-5), {{F(0.8), F(0.8), F(0.8)}, 0}},   // Back wall
    };

    // Camera
    Ray cam = {{F(0), F(0), F(0)}, {F(0), F(0), F(-1)}};

    for (int y = 0; y < height; ++y) {
        fprintf(stderr, "\rRendering: %3d%%", (y * 100) / (height - 1));
        fflush(stderr);
        for (int x = 0; x < width; ++x) {
            Vec3 total_color = {F(0), F(0), F(0)};

            for (int s = 0; s < SAMPLES_PER_PIXEL; ++s) {
                // Screen coordinates
                int32_t sx = F(2.0 * (x + I(rand_fp())) / width - 1.0);
                int32_t sy = F(1.0 - 2.0 * (y + I(rand_fp())) / height);
                
                Ray r = cam;
                r.dir.x = sx;
                r.dir.y = sy;
                r.dir = vec_norm(r.dir);

                Vec3 path_attenuation = {ONE, ONE, ONE};

                for (int b = 0; b < MAX_BOUNCES; ++b) {
                    int32_t min_t = F(1e9);
                    Material mat = {{F(0), F(0), F(0)}, 0};
                    Vec3 hit_normal = {F(0), F(0), F(0)};

                    // Sphere intersection
                    for (size_t i = 0; i < sizeof(spheres) / sizeof(Sphere); ++i) {
                        int32_t t = intersect_sphere(&r, &spheres[i]);
                        if (t < min_t) {
                            min_t = t;
                            mat = spheres[i].material;
                            Vec3 hit_point = vec_add(r.orig, vec_scale(r.dir, t));
                            hit_normal = vec_norm(vec_sub(hit_point, spheres[i].center));
                        }
                    }
                    
                    // Plane intersection
                    for (size_t i = 0; i < sizeof(planes) / sizeof(Plane); ++i) {
                        int32_t t = intersect_plane(&r, &planes[i]);
                        if (t < min_t) {
                            min_t = t;
                            mat = planes[i].material;
                            hit_normal = planes[i].normal;
                        }
                    }

                    if (min_t >= F(1e8)) {
                        break; // Ray escaped
                    }

                    if (mat.is_light) {
                        total_color = vec_add(total_color, vec_mul(path_attenuation, mat.color));
                        break;
                    }

                    // Attenuate path and prepare for next bounce
                    path_attenuation = vec_mul(path_attenuation, mat.color);

                    Vec3 hit_point = vec_add(r.orig, vec_scale(r.dir, min_t));
                    
                    // New random direction
                    Vec3 random_dir = random_unit_vector();
                    Vec3 bounce_dir = vec_add(hit_normal, random_dir);
                    if (vec_len_sq(bounce_dir) == 0) {
                        bounce_dir = hit_normal;
                    }

                    r.orig = vec_add(hit_point, vec_scale(hit_normal, F(0.001)));
                    r.dir = vec_norm(bounce_dir);
                }
            }
            
            Vec3 color = vec_scale(total_color, div_fp(ONE, F(SAMPLES_PER_PIXEL)));

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