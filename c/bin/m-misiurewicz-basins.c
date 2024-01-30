#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mandelbrot-numerics.h>
#include <mandelbrot-symbolics.h>
#include <mandelbrot-graphics.h>

void trace_ray(cairo_t *cr, m_d_transform *t, double _Complex center, double radius, const char *name, const char *s)
{
  struct m_binangle bangle;
  m_binangle_init(&bangle);
  m_binangle_from_string(&bangle, s);
  mpq_t angle;
  mpq_init(angle);
  m_binangle_to_rational(angle, &bangle);
  m_binangle_clear(&bangle);
  m_d_exray_in *ray = m_d_exray_in_new(angle, 8);
  mpq_clear(angle);
  double _Complex q = 0;
  double _Complex dc = 1;
  for (int i = 0; i < 8192; ++i)
  {
    if (m_stepped != m_d_exray_in_step(ray, 16))
    {
      break;
    }
    double _Complex c = m_d_exray_in_get(ray);
    double r = cabs(c - center) / radius;
    if (r > 0.85)
    {
fprintf(stderr, "%s %.18f %.18f\n", name, creal(q), cimag(q));
      q = c;
    }
    m_d_transform_reverse(t, &c, &dc);
    if (i > 0 && 0 < creal(c) && creal(c) < 1920 && 0 < cimag(c) && cimag(c) < 1080)
      cairo_line_to(cr, creal(c), cimag(c));
    else
      cairo_move_to(cr, creal(c), cimag(c));
  }
  cairo_stroke(cr);
  dc = 1;
fprintf(stderr, "%s %.18f %.18f\n", name, creal(q), cimag(q));
  m_d_transform_reverse(t, &q, &dc);
fprintf(stderr, "%s %f %f\n", name, creal(q), cimag(q));
  cairo_move_to(cr, creal(q), cimag(q));
  cairo_show_text(cr, name);
  cairo_fill(cr);
  m_d_exray_in_delete(ray);
}

void draw_ray(cairo_t *cr, m_d_transform *t, double _Complex center, double radius, const char *name, const char *pre, const char *per)
{
  int len = 4 + strlen(pre) + strlen(per);
  char s[len];
  s[0] = 0;
  strncat(s, ".", len);
  strncat(s, pre + 2, len);
  s[strlen(s)-1] = 0;
  strncat(s, per + 1, len);
fprintf(stderr, "%s\n", s);
  trace_ray(cr, t, center, radius, name, s);
}

void draw_ray2(cairo_t *cr, m_d_transform *t, double _Complex center, double radius, const char *name, const char *pre, const char *per1, const char *per2, const char *per3)
{
  int len = 4 + strlen(pre) + strlen(per1) + strlen(per2) + strlen(per3);
  char s[len];
  s[0] = 0;
  strncat(s, ".", len);
  strncat(s, pre + 2, len);
  s[strlen(s)-1] = 0;
  strncat(s, per1 + 1, len);
  s[strlen(s)-1] = 0;
  strncat(s, per2 + 2, len);
  s[strlen(s)-1] = 0;
  strncat(s, per3 + 2, len);
fprintf(stderr, "%s\n", s);
  trace_ray(cr, t, center, radius, name, s);
}

static inline double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

extern int main(int argc, char **argv) {
  if (argc != 10) {
    fprintf(stderr, "usage: %s out.png width height creal cimag radius maxiters preperiod period\n", argv[0]);
    return 1;
  }
  int width = atoi(argv[2]);
  int height = atoi(argv[3]);
  double _Complex center = atof(argv[4]) + I * atof(argv[5]);
  double radius = atof(argv[6]);
  int maxiters = atoi(argv[7]);
  int preperiod = atoi(argv[8]);
  int period = atoi(argv[9]);
  const char *filename = argv[1];

  m_image *img = m_image_new(width, height);
  m_d_transform *transform =
    m_d_transform_rectangular(width, height, center, radius);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < height; ++j) {
    complex double *zs = malloc(maxiters * sizeof(*zs));
    int *ps = malloc(maxiters * sizeof(*ps));
    m_compute_t bias = m_unknown;
    for (int i = 0; i < width; ++i) {
      complex double c = i + I * j;
      complex double dc0 = 1;
      m_d_transform_forward(transform, &c, &dc0);
      double pixelspacing2 = cabs2(dc0);
      double pixelspacing = sqrt(pixelspacing2);
      double de = 0;
      complex double z = c;
      complex double dc = 1;
      complex double dz = 0;
      double mz2 = 1.0 / 0.0;
      int np = 0;
      if (de == 0) {
        for (int n = 1; n < maxiters; ++n) {
          double z2 = cabs2(z);
          if (z2 < mz2) {
            mz2 = z2;
            int p = n;
            if (bias == m_interior) {
              if (m_d_interior_de(&de, &dz, z, c, p, 64)) {
                de /= pixelspacing;
                bias = m_interior;
                break;
              }
            } else {
              zs[np] = z;
              ps[np] = p;
              np++;
            }
          }
          if (z2 > 65536) {
            de = sqrt(z2 / cabs2(dc * dc0)) * log(z2);
            bias = m_exterior;
            break;
          }
          dc = 2 * z * dc + 1;
          z = z * z + c;
        }
      }
      if (de == 0 && bias != m_interior) {
        for (int n = 0; n < np; ++n) {
          if (m_d_interior_de(&de, &dz, zs[n], c, ps[n], 64)) {
            de /= pixelspacing;
            bias = m_interior;
            break;
          }
        }
      }
      m_d_misiurewicz(&z, c, preperiod, period, 64);
      m_d_misiurewicz_naive(&z, z, preperiod, period, 64);
      double hue = carg(z - center) / 6.283185307179586; hue -= floor(hue);
      double sat = cabs(z - c) < radius / 16 ? 1 : 0.25;
      double val = tanh(de + 1);
      m_image_plot(img, i, j, m_pixel_hsva(hue, sat, val, 1));
    }
    free(zs);
    free(ps);
  }
  m_image_dirty(img);
  cairo_surface_t *surface = m_image_surface(img);
  cairo_t *cr = cairo_create(surface);
  char *prelo = 0, *prehi = 0, *perlo = 0, *perhi = 0;
  double _Complex nucleus;
  m_d_nucleus(&nucleus, center, preperiod, 64);
  m_d_external_angles(&prelo, &prehi, nucleus, preperiod, 0.9, 64);
  fprintf(stderr, "%.18f %.18f\n%s\n%s\n", creal(nucleus), cimag(nucleus), prelo, prehi);
  m_d_nucleus(&nucleus, center, period, 64);
  m_d_external_angles(&perlo, &perhi, nucleus, period, 0.5, 64);
  fprintf(stderr, "%.18f %.18f\n%s\n%s\n", creal(nucleus), cimag(nucleus), perlo, perhi);
  cairo_set_source_rgba(cr, 0, 0, 0, 1);
  cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 16);
  draw_ray(cr, transform, center, radius, ".q(p)", prelo, perlo);
  draw_ray(cr, transform, center, radius, ".q(P)", prelo, perhi);
  draw_ray(cr, transform, center, radius, ".Q(p)", prehi, perlo);
  draw_ray(cr, transform, center, radius, ".Q(P)", prehi, perhi);
  //cairo_set_source_rgba(cr, 1, 0, 0, 1);
  draw_ray2(cr, transform, center, radius, ".q(ppP)", prelo, perlo, perlo, perhi);
  draw_ray2(cr, transform, center, radius, ".q(pPp)", prelo, perlo, perhi, perlo);
  //draw_ray2(cr, transform, prelo, perlo, perhi, perhi);
  //cairo_set_source_rgba(cr, 1, 1, 0, 1);
  draw_ray2(cr, transform, center, radius, ".q(Ppp)", prelo, perhi, perlo, perlo);
  //draw_ray2(cr, transform, prelo, perhi, perlo, perhi);
  //draw_ray2(cr, transform, prelo, perhi, perhi, perlo);
  //cairo_set_source_rgba(cr, 0, 1, 0, 1);
  draw_ray2(cr, transform, center, radius, ".Q(ppP)", prehi, perlo, perlo, perhi);
  draw_ray2(cr, transform, center, radius, ".Q(pPp)", prehi, perlo, perhi, perlo);
  //draw_ray2(cr, transform, prehi, perlo, perhi, perhi);
  //cairo_set_source_rgba(cr, 0, 0, 1, 1);
  draw_ray2(cr, transform, center, radius, ".Q(Ppp)", prehi, perhi, perlo, perlo);
  //draw_ray2(cr, transform, prehi, perhi, perlo, perhi);
  //draw_ray2(cr, transform, prehi, perhi, perhi, perlo);
  m_image_save_png(img, filename);
  m_d_transform_delete(transform);
  m_image_delete(img);
  return 0;
}
