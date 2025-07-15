// A C port of the GLSL shader: https://www.shadertoy.com/view/3dfGR2
// Original by Reinder Nijhoff 2019
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// Translated to C by an AI assistant.

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>

// --- Constants ---
#define PI 3.14159265359f
#define MAX_FLOAT 1e5f
#define EPSILON 0.1f
#define SAMPLES 32
#define MAX_BOUNCES 3

// Material types
#define LAMBERTIAN 0
#define DIFFUSE_LIGHT 1

// --- Vector Types and Operations ---

typedef struct { float x, y; } vec2;
typedef struct { float x, y, z; } vec3;
typedef struct { unsigned int x, y; } uvec2;
typedef struct { unsigned int x, y, z; } uvec3;

vec2 vec2_scale(vec2 a, float s) { return (vec2){a.x * s, a.y * s}; }
vec2 vec2_mul(vec2 a, vec2 b) { return (vec2){a.x * b.x, a.y * b.y}; }

vec3 vec3_add(vec3 a, vec3 b) { return (vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
vec3 vec3_sub(vec3 a, vec3 b) { return (vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
vec3 vec3_mul(vec3 a, vec3 b) { return (vec3){a.x * b.x, a.y * b.y, a.z * b.z}; }
vec3 vec3_scale(vec3 a, float s) { return (vec3){a.x * s, a.y * s, a.z * s}; }
float vec3_dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
vec3 vec3_cross(vec3 a, vec3 b) { return (vec3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
float vec3_length_sq(vec3 a) { return vec3_dot(a, a); }

// An approximation of sqrtf(x) based on a low-error piecewise linear approximation
// of the sqrt function between powers of two.
// It implements the formula: G(x) = c1*x*2^(-l/2) + c2*2^(l/2)
// where l=floor(log2(x)), c1=sqrt(2)-1, c2=2-sqrt(2).
float sqrtf_approx(float z) {
    if (z == 0.f) return 0.f;

    union { float f; int32_t i; } val = {z};
    
    // l(z) = floor(log2(z)) is approximated by extracting the exponent from the float representation.
    int32_t lz = ((val.i >> 23) & 0xFF) - 127;

    const float c1 = 0.41421356237f; // sqrt(2.f) - 1.f
    const float c2 = 0.58578643762f; // 2.f - sqrt(2.f)

    // Compute 2^(lz/2) and 2^(-lz/2) using a bit-level sqrt approximation on 2^lz and 2^-lz.
    union { int32_t i; float f; } u;

    // Compute 2^(lz/2)
    int32_t pow2_lz_i = (lz + 127) << 23;
    u.i = (pow2_lz_i >> 1) + 0x1fbd1df5;
    float pow2_half_lz = u.f;

    // Compute 2^(-lz/2)
    int32_t pow2_minus_lz_i = (-lz + 127) << 23;
    u.i = (pow2_minus_lz_i >> 1) + 0x1fbd1df5;
    float pow2_minus_half_lz = u.f;
    
    // The approximation is always an underestimation.
    return c1 * z * pow2_minus_half_lz + c2 * pow2_half_lz;
}

// An accurate sqrtf approximation that avoids division. It works by first
// computing an accurate inverse square root using 3 iterations of Newton's
// method, and then multiplying by the original number.
float sqrtf_accurate_approx(float z) {
    if (z == 0.f) return 0.f;

    int32_t i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = z * 0.5F;
    y = z;
    i = *(int32_t *)&y;
    i = 0x5f3759df - (i >> 1); // Initial guess for inverse sqrt
    y = *(float *)&i;

    // Refine the inverse sqrt guess with three iterations of Newton's method.
    y = y * (threehalfs - (x2 * y * y));
    y = y * (threehalfs - (x2 * y * y));
    y = y * (threehalfs - (x2 * y * y));

    // Multiply the inverse sqrt by the original number to get the sqrt.
    return z * y;
}

// Fast inverse square root - from Quake III Arena
// Uses a 32-bit integer to alias the float bits.
float fast_inv_sqrt(float number) {
    int32_t i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = *(int32_t *)&y;                // Evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1);          // What the f***?
    y = *(float *)&i;
    y = y * (threehalfs - (x2 * y * y)); // 1st iteration of Newton's method

    return y;
}

// Compute 1/x without division using the fast inverse square root trick
// followed by one Newton–Raphson refinement. Handles the sign of `x`.
float fast_inv(float x) {
    float absx = fabsf(x);
    if (absx == 0.f) return copysignf(MAX_FLOAT, x); // Prevent division-by-zero blow-ups
    float inv = fast_inv_sqrt(absx);
    inv *= inv;                         // ≈ 1/absx
    inv = inv * (2.f - absx * inv);     // One Newton iteration => higher accuracy
    return copysignf(inv, x);
}

vec3 vec3_normalize(vec3 a) {
    float len_sq = vec3_length_sq(a);
    if (len_sq > 0.f) {
        return vec3_scale(a, fast_inv_sqrt(len_sq));
    }
    return a;
}

// --- Random Number Generation ---

float g_seed = 0.f;
const float INV_MAX_RAND = 1.f / (float)0x7fffffffU;

unsigned int floatBitsToUint(float f) {
    union { float f; unsigned int u; } u;
    u.f = f;
    return u.u;
}

uint32_t base_hash(uvec2 p) {
    p.x = 1103515245U * ((p.x >> 1U) ^ p.y);
    p.y = 1103515245U * ((p.y >> 1U) ^ p.x);
    uint32_t h32 = 1103515245U * (p.x ^ (p.y >> 3U));
    return h32 ^ (h32 >> 16);
}

vec2 hash2(float* seed) {
    uvec2 p;
    *seed += 0.1f; p.x = floatBitsToUint(*seed);
    *seed += 0.1f; p.y = floatBitsToUint(*seed);
    uint32_t n = base_hash(p);
    uvec2 rz = { n, n * 48271U };
    return (vec2){(float)(rz.x & 0x7fffffffU) * INV_MAX_RAND, (float)(rz.y & 0x7fffffffU) * INV_MAX_RAND};
}

vec3 hash3(float* seed) {
    uvec2 p;
    *seed += 0.1f; p.x = floatBitsToUint(*seed);
    *seed += 0.1f; p.y = floatBitsToUint(*seed);
    uint32_t n = base_hash(p);
    uvec3 rz = { n, n * 16807U, n * 48271U };
    return (vec3){(float)(rz.x & 0x7fffffffU) * INV_MAX_RAND, (float)(rz.y & 0x7fffffffU) * INV_MAX_RAND, (float)(rz.z & 0x7fffffffU) * INV_MAX_RAND};
}

vec3 random_cos_weighted_hemisphere_direction(vec3 n, float* seed) {
    vec2 r = hash2(seed);
    vec3 temp = fabsf(n.y) > 0.5f ? (vec3){1.f, 0.f, 0.f} : (vec3){0.f, 1.f, 0.f};
    vec3 uu = vec3_normalize(vec3_cross(n, temp));
    vec3 vv = vec3_cross(uu, n);
    float ra = sqrtf_approx(r.y);
    float rx = ra * cosf(2.f * PI * r.x);
    float ry = ra * sinf(2.f * PI * r.x);
    float rz = sqrtf_approx(1.f - r.y);
    vec3 rr = vec3_add(vec3_add(vec3_scale(uu, rx), vec3_scale(vv, ry)), vec3_scale(n, rz));
    return vec3_normalize(rr);
}

vec2 random_in_unit_disk(float* seed) {
    vec2 h = hash2(seed);
    h.y *= 2.f * PI;
    float r = sqrtf_approx(h.x);
    return (vec2){ r * sinf(h.y), r * cosf(h.y) };
}

// --- Ray Tracing Structures ---

typedef struct { vec3 origin, direction; } ray;
typedef struct { int type; vec3 color; } material;
typedef struct { float t; vec3 p, normal; material mat; } hit_record;

// Primitives
typedef struct { vec3 center; float radius; } sphere_primitive;
typedef struct { vec3 center, dimension; } box_primitive;

// --- Material & Intersection ---

void material_scatter(ray r_in, hit_record rec, vec3* attenuation, ray* scattered) {
    scattered->origin = rec.p;
    scattered->direction = random_cos_weighted_hemisphere_direction(rec.normal, &g_seed);
    *attenuation = rec.mat.color;
}

vec3 material_emitted(hit_record rec) {
    return (rec.mat.type == DIFFUSE_LIGHT) ? rec.mat.color : (vec3){0, 0, 0};
}

// Intersectors
bool box_intersect(ray r, float t_min, float t_max, vec3 center, vec3 rad, vec3* normal, float* dist) {
    vec3 m = {fast_inv(r.direction.x), fast_inv(r.direction.y), fast_inv(r.direction.z)};
    vec3 n = vec3_mul(m, vec3_sub(r.origin, center));
    vec3 k = vec3_mul((vec3){fabsf(m.x), fabsf(m.y), fabsf(m.z)}, rad);
    vec3 t1 = vec3_sub((vec3){-n.x, -n.y, -n.z}, k);
    vec3 t2 = vec3_add((vec3){-n.x, -n.y, -n.z}, k);
    float tN = fmaxf(fmaxf(t1.x, t1.y), t1.z);
    float tF = fminf(fminf(t2.x, t2.y), t2.z);
    if (tN > tF || tF < 0.f) return false;
    float t = tN < t_min ? tF : tN;
    if (t < t_max && t > t_min) {
        *dist = t;
        if (t1.x > t1.y && t1.x > t1.z) *normal = (vec3){ r.direction.x > 0.f ? -1.f : 1.f, 0, 0};
        else if (t1.y > t1.z) *normal = (vec3){0, r.direction.y > 0.f ? -1.f : 1.f, 0};
        else *normal = (vec3){0, 0, r.direction.z > 0.f ? -1.f : 1.f};
        return true;
    }
    return false;
}

bool sphere_hit(sphere_primitive s, ray r, float t_min, float t_max, hit_record* rec) {
    vec3 oc = vec3_sub(r.origin, s.center);
    float a = vec3_dot(r.direction, r.direction);
    float inv_a = fast_inv(a);
    float half_b = vec3_dot(oc, r.direction);
    float c = vec3_dot(oc, oc) - s.radius*s.radius;
    float discriminant = half_b*half_b - a*c;

    if (discriminant < 0) return false;
    
    float sqrtd = sqrtf_accurate_approx(discriminant);
    float root = (-half_b - sqrtd) * inv_a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) * inv_a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec->t = root;
    rec->p = vec3_add(r.origin, vec3_scale(r.direction, rec->t));
    vec3 outward_normal = vec3_scale(vec3_sub(rec->p, s.center), fast_inv(s.radius));
    rec->normal = outward_normal;

    return true;
}

// Shadow-specific intersection tests that check against squared distance
bool box_intersect_shadow(ray r, float t_min, float t_max_sq, vec3 center, vec3 rad) {
    vec3 m = {fast_inv(r.direction.x), fast_inv(r.direction.y), fast_inv(r.direction.z)};
    vec3 n = vec3_mul(m, vec3_sub(r.origin, center));
    vec3 k = vec3_mul((vec3){fabsf(m.x), fabsf(m.y), fabsf(m.z)}, rad);
    vec3 t1 = vec3_sub((vec3){-n.x, -n.y, -n.z}, k);
    vec3 t2 = vec3_add((vec3){-n.x, -n.y, -n.z}, k);
    float tN = fmaxf(fmaxf(t1.x, t1.y), t1.z);
    float tF = fminf(fminf(t2.x, t2.y), t2.z);
    if (tN > tF || tF < 0.f) return false;
    float t = tN < t_min ? tF : tN;
    if (t > t_min && t*t < t_max_sq) {
        return true;
    }
    return false;
}

bool sphere_hit_shadow(sphere_primitive s, ray r, float t_min, float t_max_sq) {
    vec3 oc = vec3_sub(r.origin, s.center);
    float a = vec3_dot(r.direction, r.direction);
    float inv_a = fast_inv(a);
    float half_b = vec3_dot(oc, r.direction);
    float c = vec3_dot(oc, oc) - s.radius*s.radius;
    float discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    float sqrtd = sqrtf_accurate_approx(discriminant);
    float root = (-half_b - sqrtd) * inv_a;
    if (root > t_min && root*root < t_max_sq) {
        return true;
    }
    root = (-half_b + sqrtd) * inv_a;
    if (root > t_min && root*root < t_max_sq) {
        return true;
    }
    return false;
}


// --- Camera ---

typedef struct { vec3 origin, lower_left_corner, horizontal, vertical, u, v, w; float lens_radius; } camera;

camera camera_const(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect, float aperture, float focus_dist) {
    camera cam;
    cam.lens_radius = aperture / 2.f;
    float theta = vfov * PI / 180.f;
    float half_height = tanf(theta / 2.f);
    float half_width = aspect * half_height;
    cam.origin = lookfrom;
    cam.w = vec3_normalize(vec3_sub(lookfrom, lookat));
    cam.u = vec3_normalize(vec3_cross(vup, cam.w));
    cam.v = vec3_cross(cam.w, cam.u);
    cam.lower_left_corner = vec3_sub(vec3_sub(vec3_sub(cam.origin, vec3_scale(cam.u, half_width * focus_dist)), vec3_scale(cam.v, half_height * focus_dist)), vec3_scale(cam.w, focus_dist));
    cam.horizontal = vec3_scale(cam.u, 2.f * half_width * focus_dist);
    cam.vertical = vec3_scale(cam.v, 2.f * half_height * focus_dist);
    return cam;
}
    
ray camera_get_ray(camera c, vec2 uv) {
    vec2 rd = vec2_scale(random_in_unit_disk(&g_seed), c.lens_radius);
    vec3 offset = vec3_add(vec3_scale(c.u, rd.x), vec3_scale(c.v, rd.y));
    vec3 origin = vec3_add(c.origin, offset);
    vec3 target = vec3_add(vec3_add(c.lower_left_corner, vec3_scale(c.horizontal, uv.x)), vec3_scale(c.vertical, uv.y));
    return (ray){origin, vec3_normalize(vec3_sub(target, origin))};
}

// --- Scene Definition ---

material mat_red, mat_white, mat_green, mat_light, mat_sphere1, mat_sphere2;
box_primitive box_floor, box_ceiling, box_back_wall, box_right_wall, box_left_wall, box_light;
sphere_primitive sphere1, sphere2;

void initialize_scene() {
    mat_red = (material){LAMBERTIAN, {0.65f, 0.05f, 0.05f}};
    mat_white = (material){LAMBERTIAN, {0.73f, 0.73f, 0.73f}};
    mat_green = (material){LAMBERTIAN, {0.12f, 0.45f, 0.15f}};
    mat_light = (material){DIFFUSE_LIGHT, {15.f, 15.f, 15.f}};
    mat_sphere1 = (material){LAMBERTIAN, {0.8f, 0.8f, 0.8f}};
    mat_sphere2 = (material){LAMBERTIAN, {0.8f, 0.6f, 0.2f}};

    box_floor = (box_primitive){{277.5, 0, 277.5}, {277.5, 1, 277.5}};
    box_ceiling = (box_primitive){{277.5, 555, 277.5}, {277.5, 1, 277.5}};
    box_back_wall = (box_primitive){{277.5, 277.5, 555}, {277.5, 277.5, 1}};
    box_right_wall = (box_primitive){{0, 277.5, 277.5}, {1, 277.5, 277.5}};
    box_left_wall = (box_primitive){{555, 277.5, 277.5}, {1, 277.5, 277.5}};
    box_light = (box_primitive){{278, 554, 279.5}, {65, 1, 52.5}};

    sphere1 = (sphere_primitive){{190, 90, 190}, 90};
    sphere2 = (sphere_primitive){{400, 150, 230}, 110};
}


bool world_hit(ray r, float t_min, hit_record* rec) {
    hit_record temp_rec;
    bool hit_anything = false;
    float closest_so_far = MAX_FLOAT;
    vec3 normal;
    float dist;

    if (box_intersect(r, t_min, closest_so_far, box_floor.center, box_floor.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_white;
    }
    if (box_intersect(r, t_min, closest_so_far, box_ceiling.center, box_ceiling.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_white;
    }
    if (box_intersect(r, t_min, closest_so_far, box_back_wall.center, box_back_wall.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_white;
    }
    if (box_intersect(r, t_min, closest_so_far, box_right_wall.center, box_right_wall.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_green;
    }
    if (box_intersect(r, t_min, closest_so_far, box_left_wall.center, box_left_wall.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_red;
    }
    if (box_intersect(r, t_min, closest_so_far, box_light.center, box_light.dimension, &normal, &dist)) {
        hit_anything = true; closest_so_far = dist; rec->t = dist; rec->normal = normal; rec->mat = mat_light;
    }
    if (sphere_hit(sphere1, r, t_min, closest_so_far, &temp_rec)) {
        hit_anything = true; closest_so_far = temp_rec.t; rec->t = temp_rec.t; rec->normal = temp_rec.normal; rec->mat = mat_sphere1;
    }
    if (sphere_hit(sphere2, r, t_min, closest_so_far, &temp_rec)) {
        hit_anything = true; closest_so_far = temp_rec.t; rec->t = temp_rec.t; rec->normal = temp_rec.normal; rec->mat = mat_sphere2;
    }

    if(hit_anything) {
        rec->p = vec3_add(r.origin, vec3_scale(r.direction, rec->t));
    }
    return hit_anything;
}

bool shadow_hit(ray r, float t_min, float t_max_sq) {
    if (box_intersect_shadow(r, t_min, t_max_sq, box_floor.center, box_floor.dimension)) return true;
    if (box_intersect_shadow(r, t_min, t_max_sq, box_ceiling.center, box_ceiling.dimension)) return true;
    if (box_intersect_shadow(r, t_min, t_max_sq, box_back_wall.center, box_back_wall.dimension)) return true;
    if (box_intersect_shadow(r, t_min, t_max_sq, box_right_wall.center, box_right_wall.dimension)) return true;
    if (box_intersect_shadow(r, t_min, t_max_sq, box_left_wall.center, box_left_wall.dimension)) return true;
    if (sphere_hit_shadow(sphere1, r, t_min, t_max_sq)) return true;
    if (sphere_hit_shadow(sphere2, r, t_min, t_max_sq)) return true;
    
    return false;
}

vec3 color(ray r) {
    vec3 col = {1.f, 1.f, 1.f}; // Throughput
    vec3 emitted = {0.f, 0.f, 0.f};
    hit_record rec;

    for (int bounce = 0; bounce < MAX_BOUNCES; bounce++) {
        rec.t = MAX_FLOAT;
        if (!world_hit(r, EPSILON, &rec)) {
            return emitted;
        }
        if (rec.mat.type == DIFFUSE_LIGHT) {
            if (bounce == 0) {
                return rec.mat.color;
            } else {
                return emitted;
            }
        }

        vec3 attenuation;
        ray scattered;
        material_scatter(r, rec, &attenuation, &scattered);
        if (bounce == 0) {
            col = attenuation;
        } else {
            col = vec3_mul(col, attenuation);
        }
        r = scattered;

        // Direct light sampling for bounce
        vec3 p_on_light = vec3_add(box_light.center,
            vec3_mul(vec3_sub(vec3_scale(hash3(&g_seed), 2.f), (vec3){1,1,1}), box_light.dimension));
        vec3 L = vec3_sub(p_on_light, rec.p);
        float rr = vec3_dot(L, L);
        L = vec3_normalize(L);

        vec3 origin_offset = vec3_add(rec.p, vec3_scale(rec.normal, 1));
        ray shadowRay = {origin_offset, L};
        float rr_shadow = vec3_length_sq(vec3_sub(p_on_light, origin_offset));

        if (vec3_dot(rec.normal, L) > 0.f && L.y > 0.f && !shadow_hit(shadowRay, .01, rr_shadow)) {
            const float light_area = (box_light.dimension.x * box_light.dimension.z) * 4.f;
            float denom = (PI * rr + EPSILON);
            float inv_denom = fast_inv_sqrt(denom);
            inv_denom *= inv_denom;
            float weight = light_area * L.y * vec3_dot(rec.normal, L) * inv_denom;
            emitted = vec3_add(emitted,
                vec3_mul(col, vec3_scale((vec3){15.f, 15.f, 15.f}, weight)));
        }
    }

    return emitted;
}


// --- Main ---

int main() {
    initialize_scene();
    int width = 512;
    int height = 512;
    float aspect = (float)width / (float)height;

    vec3 lookfrom = {278.f, 278.f, -800.f};
    vec3 lookat = {278, 278, 0};

    FILE* f = fopen("image.ppm", "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "P3\n%d %d\n255\n", width, height);

    for (int j = height - 1; j >= 0; j--) {
        fprintf(stderr, "\rScanlines remaining: %d ", j);
        for (int i = 0; i < width; i++) {
            
            vec3 tcol = {0,0,0};
            
            for (int s=0; s < SAMPLES; s++) {
                g_seed = (float)j*width + i + (float)s/SAMPLES;
                
                vec2 rand_offset = hash2(&g_seed);
                vec2 uv = {((float)i + rand_offset.x) / (float)width, ((float)j + rand_offset.y) / (float)height};
                
                camera cam = camera_const(lookfrom, lookat, (vec3){0,1,0}, 40.f, aspect, .0f, 10.f);
                ray r = camera_get_ray(cam, uv);
                tcol = vec3_add(tcol, color(r));
            }
            
            tcol = vec3_scale(tcol, 1.f / (float)SAMPLES);
            tcol.x = sqrtf(tcol.x); tcol.y = sqrtf(tcol.y); tcol.z = sqrtf(tcol.z);
            
            int ir = fminf(255.f, 255.99f * tcol.x);
            int ig = fminf(255.f, 255.99f * tcol.y);
            int ib = fminf(255.f, 255.99f * tcol.z);
            fprintf(f, "%d %d %d\n", ir, ig, ib);
        }
    }
    fprintf(stderr, "\nDone.\n");
    fclose(f);
    return 0;
} 