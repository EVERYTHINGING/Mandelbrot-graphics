#include <mandelbrot-numerics.h>
#include <mandelbrot-graphics.h>
#include <stdio.h>

static const double twopi = 6.283185307179586;

extern int main(int argc, char **argv) {
  int w = 256 * 6;
  int h = 256;
  double r = 2.0/3.0;
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 100;
  int maxiters = 8192;
  const char *filename = "out.png";

  if (! (argc == 7)) {
    return 1;
  }
  complex double nucleus = atof(argv[1]) + I * atof(argv[2]);
  int period = atoi(argv[3]);
  mpq_t tzero, tone, tinfinity;
  mpq_init(tzero);
  mpq_init(tone);
  mpq_init(tinfinity);
  mpq_set_str(tzero, argv[4], 10);
  mpq_set_str(tone, argv[5], 10);
  mpq_set_str(tinfinity, argv[6], 10);
  mpq_canonicalize(tzero);
  mpq_canonicalize(tone);
  mpq_canonicalize(tinfinity);

  m_d_nucleus(&nucleus, nucleus, period, 64);
  complex double zero, one, infinity, z;
  m_d_interior(&z, &zero, nucleus, nucleus, cexp(I * twopi * mpq_get_d(tzero)), period, 64);
  m_d_interior(&z, &one, nucleus, nucleus, cexp(I * twopi * mpq_get_d(tone)), period, 64);
  if (mpq_sgn(tinfinity) == 0) {
    m_d_parent(tinfinity, &infinity, &z, nucleus, period, 64);
  } else {
    m_d_interior(&z, &infinity, nucleus, nucleus, cexp(I * twopi * mpq_get_d(tinfinity)), period, 64);
  }

  printf("%.16e %.16e\n", creal(zero), cimag(zero));
  printf("%.16e %.16e\n", creal(one), cimag(one));
  printf("%.16e %.16e\n", creal(infinity), cimag(infinity));

  m_d_transform *rect = m_d_transform_rectangular(w, h, I * 0.99 * r, r);
  m_d_transform *transform = 0;
  switch (m_d_shape(nucleus, period)) {
    case m_cardioid: {
      complex double size = m_d_size(nucleus, period);
      complex double half, cusp;
      m_d_interior(&z, &half, nucleus, nucleus, -1, period, 64);
      m_d_interior(&z, &cusp, nucleus, nucleus, 1, period, 64);
//      complex double size = cusp - half;
      m_d_transform *cardioid = m_d_transform_cardioid();
      m_d_transform *linear = m_d_transform_linear(-nucleus / size, 1 / size);
      m_d_transform *t0 = m_d_transform_compose(cardioid, linear);
      m_d_transform_reverse(t0, &zero, &z);
      m_d_transform_reverse(t0, &one, &z);
      m_d_transform_reverse(t0, &infinity, &z);
      m_d_transform *moebius = m_d_transform_moebius3(zero, one, infinity);
      m_d_transform *t1 = m_d_transform_compose(rect, moebius);
      transform = m_d_transform_compose(t1, t0);
      break;
    }
    case m_circle: {
      m_d_transform *moebius = m_d_transform_moebius3(zero, one, infinity);
      transform = m_d_transform_compose(rect, moebius);
      break;
    }
    default: {
      return 1;
    }
  }

  m_d_colour_t *colour = m_d_colour_minimal(white, black, white);
  m_image *image = m_image_new(w, h);
  m_d_render_scanline(image, transform, er, maxiters, colour);
  m_image_save_png(image, filename);

  return 0;
}
