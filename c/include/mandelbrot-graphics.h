#ifndef MANDELBROT_GRAPHICS_H
#define MANDELBROT_GRAPHICS_H 1

#include <cairo.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include <mandelbrot-numerics.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* functions returning bool return true for success, false for failure */

enum m_compute_t {
  m_unknown,
  m_interior,
  m_exterior
};
typedef enum m_compute_t m_compute_t;

/* double precision: m_d_*()  */

struct m_d_compute;
typedef struct m_d_compute m_d_compute;

extern m_d_compute *m_d_compute_alloc(int npartials);
extern void m_d_compute_free(m_d_compute *px);
extern void m_d_compute_init(m_d_compute *px, m_compute_t bias, double er, double _Complex c, m_period_filter_t *filter);
extern void m_d_compute_clear(m_d_compute *px);
extern bool m_d_compute_step(m_d_compute *px, int steps);
extern m_compute_t m_d_compute_get_tag(const m_d_compute *px);
extern int m_d_compute_get_n(const m_d_compute *px);
extern int m_d_compute_get_p(const m_d_compute *px);
extern int m_d_compute_get_q(const m_d_compute *px);
extern double _Complex m_d_compute_get_c(const m_d_compute *px);
extern double _Complex m_d_compute_get_z(const m_d_compute *px);
extern double _Complex m_d_compute_get_dc(const m_d_compute *px);
extern double _Complex m_d_compute_get_dz(const m_d_compute *px);
extern double _Complex m_d_compute_get_zp(const m_d_compute *px);
extern double _Complex m_d_compute_get_zq(const m_d_compute *px);
extern double m_d_compute_get_de(const m_d_compute *px);

struct m_d_transform;
typedef struct m_d_transform m_d_transform;
typedef void (m_d_transform_f)(void *userdata, double _Complex *c, double _Complex *dc);
struct m_d_transform {
  m_d_transform_f *forward;
  m_d_transform_f *reverse;
};

extern void m_d_transform_forward(m_d_transform *t, double _Complex *c, double _Complex *dc);
extern void m_d_transform_reverse(m_d_transform *t, double _Complex *c, double _Complex *dc);
extern void m_d_transform_delete(m_d_transform *t);
extern m_d_transform *m_d_transform_compose(m_d_transform *s, m_d_transform *t);
extern m_d_transform *m_d_transform_invert(m_d_transform *t);
extern m_d_transform *m_d_transform_rectangular(int w, int h, double _Complex c, double r);
extern m_d_transform *m_d_transform_linear(double _Complex add, double _Complex mul);
extern m_d_transform *m_d_transform_moebius(const m_d_mat2 *m);
extern m_d_transform *m_d_transform_moebius3(double _Complex zero, double _Complex one, double _Complex infinity);
extern m_d_transform *m_d_transform_moebius6(double _Complex a0, double _Complex a1, double _Complex a2, double _Complex b0, double _Complex b1, double _Complex b2);
extern m_d_transform *m_d_transform_cardioid();
extern m_d_transform *m_d_transform_exponential(double _Complex center);

typedef uint32_t m_pixel_t;

extern double m_pixel_red(m_pixel_t p);
extern double m_pixel_green(m_pixel_t p);
extern double m_pixel_blue(m_pixel_t p);
extern double m_pixel_alpha(m_pixel_t p);
extern m_pixel_t m_pixel_rgba(double r, double g, double b, double a);
extern m_pixel_t m_pixel_hsva(double h, double s, double v, double a);
extern m_pixel_t m_pixel_mix(m_pixel_t a, m_pixel_t b, double x);
extern m_pixel_t m_pixel_blacken(m_pixel_t a, double x);
extern m_pixel_t m_pixel_whiten(m_pixel_t a, double x);

struct m_image;
typedef struct m_image m_image;

extern m_image *m_image_new(int width, int height);
extern void m_image_delete(m_image *img);
extern m_image *m_image_load_ppm(const char *filename);
extern void m_image_flush(m_image *img);
extern void m_image_dirty(m_image *img);
extern bool m_image_in_bounds(m_image *img, int x, int y);
extern void m_image_plot(m_image *img, int x, int y, m_pixel_t c);
extern m_pixel_t m_image_peek(m_image *img, int x, int y);
extern bool m_image_save_png(m_image *img, const char *filename);
extern int m_image_get_width(const m_image *img);
extern int m_image_get_height(const m_image *img);
extern cairo_surface_t *m_image_surface(m_image *img);

struct m_mipmap;
typedef struct m_mipmap m_mipmap;

extern m_mipmap *m_mipmap_new(m_image *img);
extern m_pixel_t m_mipmap_linear(const m_mipmap *m, double level, double x, double y);
extern m_pixel_t m_mipmap_linear_linear(const m_mipmap *m, double level, double x, double y);
extern int m_mipmap_get_levels(const m_mipmap *m);

struct m_d_colour_t;
typedef struct m_d_colour_t m_d_colour_t;
typedef m_pixel_t (m_d_colour_f)(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc);
struct m_d_colour_t {
  m_d_colour_f *colour;
};

extern void m_d_colour_delete(m_d_colour_t *colour);
extern m_pixel_t m_d_colour(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc);
extern m_d_colour_t *m_d_colour_minimal(m_pixel_t interior, m_pixel_t boundary, m_pixel_t exterior);
extern double _Complex m_d_colour_exterior_coordinates(m_d_compute *px, double er);
extern double m_d_colour_grid(m_d_compute *px, double er);
extern m_d_colour_t *m_d_colour_minimal_grid(m_pixel_t interior, m_pixel_t boundary, m_pixel_t exterior, m_pixel_t grid);
extern m_d_colour_t *m_d_colour_minimal_texture(m_pixel_t interior, m_pixel_t boundary, m_mipmap *exterior);
extern m_d_colour_t *m_d_colour_domain(double hue, double sat, double val, m_period_filter_t *filter);

extern void m_d_render_scanline_filtered(m_image *image, m_d_transform *transform, double er, int maxiters, m_d_colour_t *colour, m_period_filter_t *filter);
extern void m_d_render_scanline(m_image *image, m_d_transform *transform, double er, int maxiters, m_d_colour_t *colour);

#ifdef __cplusplus
}
#endif

#endif
