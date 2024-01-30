#include <mandelbrot-graphics.h>
#include "m_d_util.h"

extern double m_pixel_red(m_pixel_t p) {
  return ((p & 0x00ff0000) >> 16) / (double) ((p & 0xff000000) >> 24);
}

extern double m_pixel_green(m_pixel_t p) {
  return ((p & 0x0000ff00) >> 8) / (double) ((p & 0xff000000) >> 24);
}

extern double m_pixel_blue(m_pixel_t p) {
  return ((p & 0x000000ff)) / (double) ((p & 0xff000000) >> 24);
}

extern double m_pixel_alpha(m_pixel_t p) {
  return ((p & 0xff000000) >> 24) / 255.0;
}

extern m_pixel_t m_pixel_rgba(double r, double g, double b, double a) {
  m_pixel_t ri = fmin(fmax(255 * r * a + 0.5, 0), 255);
  m_pixel_t gi = fmin(fmax(255 * g * a + 0.5, 0), 255);
  m_pixel_t bi = fmin(fmax(255 * b * a + 0.5, 0), 255);
  m_pixel_t ai = fmin(fmax(255 * a + 0.5, 0), 255);
  return (ai << 24) | (ri << 16) | (gi << 8) | bi;
}

extern m_pixel_t m_pixel_hsva(double h, double s, double v, double a) {
  double i, f, p, q, t, r, g, b;
  int ii;
  if (s == 0.0) { r = g = b = v; } else {
    h = 6 * (h - floor(h));
    ii = i = floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
      case 0: r = v; g = t; b = p; break;
      case 1: r = q; g = v; b = p; break;
      case 2: r = p; g = v; b = t; break;
      case 3: r = p; g = q; b = v; break;
      case 4: r = t; g = p; b = v; break;
      default:r = v; g = p; b = q; break;
    }
  }
  return m_pixel_rgba(r, g, b, a);
}

extern m_pixel_t m_pixel_mix(m_pixel_t a, m_pixel_t b, double x) {
  return m_pixel_rgba
    ( mix(m_pixel_red(a), m_pixel_red(b), x)
    , mix(m_pixel_green(a), m_pixel_green(b), x)
    , mix(m_pixel_blue(a), m_pixel_blue(b), x)
    , mix(m_pixel_alpha(a), m_pixel_alpha(b), x)
    );
}

extern m_pixel_t m_pixel_blacken(m_pixel_t a, double x) {
  return m_pixel_mix(a, m_pixel_rgba(0, 0, 0, 1), x);
}

extern m_pixel_t m_pixel_whiten(m_pixel_t a, double x) {
  return m_pixel_mix(a, m_pixel_rgba(1, 1, 1, 1), x);
}
