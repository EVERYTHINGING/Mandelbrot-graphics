#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include <mandelbrot-numerics.h>
#include <mandelbrot-graphics.h>

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s output.png angle...\n"
    , progname
    );
}

static bool arg_rational(const char *arg, mpq_t x) {
  int ok = mpq_set_str(x, arg, 10);
  mpq_canonicalize(x);
  return ok == 0;
}

static const double twopi = 6.283185307179586;

extern int main(int argc, char **argv) {
  if (! (argc > 2)) {
    usage(argv[0]);
    return 1;
  }
  const int maxsteps = 64;
  const int sharpness = 256;
  mpq_t q;
  mpq_init(q);
  double dperiod = 1;
  for (int arg = 2; arg < argc; ++arg)
  {
    if (! arg_rational(argv[arg], q)) { mpq_clear(q); return 1; }
    if (! (0 < mpq_cmp_si(q, 0, 1) && mpq_cmp_si(q, 1, 1) < 1)) { mpq_clear(q); return 1; }
    dperiod *= mpz_get_d(mpq_denref(q));
  }
  if (dperiod > INT_MAX) { mpq_clear(q); return 1; }
  int height = 1024;
  int width  = height;
  m_image *img = m_image_new(width, height);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  for (int y = 0; y < height; ++y)
  for (int x = 0; x < width; ++x)
    m_image_plot(img, x, y, black);
  int period = 1;
  double _Complex rootc = 0.25;
  double _Complex rootz = 0.5;
  double _Complex rooti = 1;
  for (int arg = 2; arg < argc; ++arg)
  {
    arg_rational(argv[arg], q);
    double t = twopi * mpq_get_d(q);

    double _Complex nucleus;
    m_d_nucleus(&nucleus, rootc, period, maxsteps);
    double _Complex bondi = cexp(I * t);
    double _Complex bondc = 0;
    double _Complex bondz = 0;
    m_d_interior(&bondz, &bondc, 0, nucleus, bondi, period, maxsteps);
    for (int step = 0; step < sharpness; ++step)
    {
      double k = (step + 0.5) / (sharpness);
      double k1 = 1 - k;
      double _Complex c = rootc * k1 + k * bondc;
      double _Complex z = rootz * k1 + k * bondz;
      double _Complex i = rooti * k1 + k * bondi;
      m_d_interior(&z, &c, 0, nucleus, i, period, maxsteps);
      double _Complex w = z;
      for (int att = 0; att < period; ++att)
      {
        m_pixel_t colour = m_pixel_hsva(att / (double) period, 1, 1, 1);
        int x = (2 + creal(w)) * width / 4;
        int y = (2 - cimag(w)) * height / 4;
        m_image_plot(img, x, y, colour);
        w = w * w + c;
      }
    }
    period *= mpz_get_ui(mpq_denref(q));;
    rootc = bondc;
    rootz = bondz;
  }
  m_image_dirty(img);
  m_image_save_png(img, argv[1]);
  m_image_delete(img);
  mpq_clear(q);
  return 0;
}
