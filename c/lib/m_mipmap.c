#include <mandelbrot-graphics.h>
#include "m_d_util.h"

struct m_mipmap {
  int levels;
  m_image **images;
};

extern m_mipmap *m_mipmap_new(m_image *img) {
  int levels = 1 + round(log2(m_image_get_width(img)));
  int size = 1 << (levels - 1);
  if (! (m_image_get_width(img) == size && m_image_get_height(img) == size)) {
    return 0;
  }
  m_mipmap *m = malloc(sizeof(*m));
  m->levels = levels;
  m->images = malloc(m->levels * sizeof(*m->images));
  m->images[0] = img;
  for (int l = 1; l < m->levels; ++l) {
    int w = m_image_get_width(m->images[l-1])/2;
    int h = m_image_get_height(m->images[l-1])/2;
    m->images[l] = m_image_new(w, h);
    #pragma omp parallel for
    for (int j = 0; j < h; ++j) {
      for (int i = 0; i < w; ++i) {
        m_pixel_t p00 = m_image_peek(m->images[l-1], 2 * i + 0, 2 * j + 0);
        m_pixel_t p10 = m_image_peek(m->images[l-1], 2 * i + 1, 2 * j + 0);
        m_pixel_t p01 = m_image_peek(m->images[l-1], 2 * i + 0, 2 * j + 1);
        m_pixel_t p11 = m_image_peek(m->images[l-1], 2 * i + 1, 2 * j + 1);
        m_pixel_t p0 = m_pixel_mix(p00, p01, 0.5);
        m_pixel_t p1 = m_pixel_mix(p10, p11, 0.5);
        m_pixel_t p = m_pixel_mix(p0, p1, 0.5);
        m_image_plot(m->images[l], i, j, p);
      }
    }
    m_image_dirty(m->images[l]);
  }
  return m;
}

extern m_pixel_t m_mipmap_linear(const m_mipmap *m, double level, double x, double y) {
  int l = clamp(round(level), 0, m->levels - 1);
  int w = m_image_get_width(m->images[l]);
  int h = m_image_get_width(m->images[l]);
  double i = fmod(fmod(w * x, w) + w, w);
  double j = h * y;
  int i0 = floor(i);
  int j0 = floor(j);
  int i1 = i0 + 1;
  int j1 = j0 + 1;
  double ix = i - i0;
  double jx = j - j0;
  i0 = ((i0 % w) + w) % w;
  j0 = clamp(j0, 0, h - 1);
  i1 = ((i1 % w) + w) % w;
  j1 = clamp(j1, 0, h - 1);
  m_pixel_t p00 = m_image_peek(m->images[l], i0, j0);
  m_pixel_t p10 = m_image_peek(m->images[l], i1, j0);
  m_pixel_t p01 = m_image_peek(m->images[l], i0, j1);
  m_pixel_t p11 = m_image_peek(m->images[l], i1, j1);
  m_pixel_t p0 = m_pixel_mix(p00, p01, jx);
  m_pixel_t p1 = m_pixel_mix(p10, p11, jx);
  m_pixel_t p = m_pixel_mix(p0, p1, ix);
  return p;
}

extern m_pixel_t m_mipmap_linear_linear(const m_mipmap *m, double level, double x, double y) {
  int l0 = floor(level);
  int l1 = l0 + 1;
  double lx = level - l0;
  return m_pixel_mix(m_mipmap_linear(m, l0, x, y), m_mipmap_linear(m, l1, x, y), lx);
}

extern int m_mipmap_get_levels(const m_mipmap *m) {
  return m->levels;
}
