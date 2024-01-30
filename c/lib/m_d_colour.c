#include <mandelbrot-graphics.h>
#include "m_d_util.h"

extern void m_d_colour_delete(m_d_colour_t *colour) {
  free(colour);
}

extern m_pixel_t m_d_colour(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc) {
  return colour->colour(colour, px, er, dc);
}


struct m_d_colour_minimal_t {
  m_d_colour_t colour;
  m_pixel_t interior;
  m_pixel_t boundary;
  m_pixel_t exterior;
};
typedef struct m_d_colour_minimal_t m_d_colour_minimal_t;

static m_pixel_t m_d_colour_minimal_do(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc) {
  (void) er;
  m_d_colour_minimal_t *c = (m_d_colour_minimal_t *) colour;
  m_pixel_t base = c->boundary;
  switch (m_d_compute_get_tag(px)) {
    case m_exterior: base = c->exterior; break;
    case m_interior: base = c->interior; break;
    default:         base = c->boundary; break;
  }
  double de = clamp(m_d_compute_get_de(px) / cabs(dc), 0, 8);
  return m_pixel_mix(c->boundary, base, tanh(de));
}

extern m_d_colour_t *m_d_colour_minimal(m_pixel_t interior, m_pixel_t boundary, m_pixel_t exterior) {
  m_d_colour_minimal_t *c = malloc(sizeof(*c));
  c->colour.colour = m_d_colour_minimal_do;
  c->interior = interior;
  c->boundary = boundary;
  c->exterior = exterior;
  return &c->colour;
}


extern double _Complex m_d_colour_exterior_coordinates(m_d_compute *px, double er) {
  double _Complex z = m_d_compute_get_z(px);
  double r = log(cabs(z)) / log(er) - 1.0;
  double a = fmod(carg(z) / twopi + 1.0, 1.0);
  return a + I * r;
}


extern double m_d_colour_grid(m_d_compute *px, double er) {
  double _Complex z = m_d_colour_exterior_coordinates(px, er);
  double r = cimag(z);
  double a = creal(z);
  double k = pow(0.5, 0.5 - r);
  double w = 0.05;
  return 1.0 - (w < r && r < 1.0 - w && w * k < a && a < 1.0 - w * k);
}


struct m_d_colour_minimal_grid_t {
  m_d_colour_t colour;
  m_pixel_t interior;
  m_pixel_t boundary;
  m_pixel_t exterior;
  m_pixel_t grid;
};
typedef struct m_d_colour_minimal_grid_t m_d_colour_minimal_grid_t;

static m_pixel_t m_d_colour_minimal_grid_do(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc) {
  (void) er;
  m_d_colour_minimal_grid_t *c = (m_d_colour_minimal_grid_t *) colour;
  m_pixel_t base = c->boundary;
  double grid = 0.0;
  switch (m_d_compute_get_tag(px)) {
    case m_exterior: base = c->exterior; grid = m_d_colour_grid(px, er); break;
    case m_interior: base = c->interior; break;
    default:         base = c->boundary; break;
  }
  base = m_pixel_mix(base, c->grid, grid);
  double de = clamp(m_d_compute_get_de(px) / cabs(dc), 0, 8);
  return m_pixel_mix(c->boundary, base, tanh(de));
}

extern m_d_colour_t *m_d_colour_minimal_grid(m_pixel_t interior, m_pixel_t boundary, m_pixel_t exterior, m_pixel_t grid) {
  m_d_colour_minimal_grid_t *c = malloc(sizeof(*c));
  c->colour.colour = m_d_colour_minimal_grid_do;
  c->interior = interior;
  c->boundary = boundary;
  c->exterior = exterior;
  c->grid = grid;
  return &c->colour;
}


struct m_d_colour_minimal_texture_t {
  m_d_colour_t colour;
  m_pixel_t interior;
  m_pixel_t boundary;
  m_mipmap *exterior;
};
typedef struct m_d_colour_minimal_texture_t m_d_colour_minimal_texture_t;


static m_pixel_t m_d_colour_minimal_texture_do(m_d_colour_t *colour, m_d_compute *px, double er, double _Complex dc) {
  m_d_colour_minimal_texture_t *c = (m_d_colour_minimal_texture_t *) colour;
  double de = m_d_compute_get_de(px) / cabs(dc);
  m_pixel_t base = c->boundary;
  switch (m_d_compute_get_tag(px)) {
    case m_exterior: {
      double level = m_mipmap_get_levels(c->exterior) - log2(de);
      double _Complex z = m_d_colour_exterior_coordinates(px, er);
      base = m_mipmap_linear_linear(c->exterior, level, creal(z), 1 - cimag(z));
      break;
    }
    case m_interior: base = c->interior; break;
    default:         base = c->boundary; break;
  }
  return m_pixel_mix(c->boundary, base, tanh(clamp(de, 0, 8)));
}

extern m_d_colour_t *m_d_colour_minimal_texture(m_pixel_t interior, m_pixel_t boundary, m_mipmap *exterior) {
  m_d_colour_minimal_texture_t *c = malloc(sizeof(*c));
  c->colour.colour = m_d_colour_minimal_texture_do;
  c->interior = interior;
  c->boundary = boundary;
  c->exterior = exterior;
  return &c->colour;
}


struct m_d_colour_domain_t {
  m_d_colour_t colour;
  double hue;
  double sat;
  double val;
  m_period_filter_t *filter;
};
typedef struct m_d_colour_domain_t m_d_colour_domain_t;

static m_pixel_t m_d_colour_domain_do(m_d_colour_t *c0, m_d_compute *px, double er, double _Complex dc) {
  (void) er;
  m_d_colour_domain_t *c = (m_d_colour_domain_t *) c0;
  int p = m_d_compute_get_p(px);
  double _Complex zp = m_d_compute_get_zp(px);
  if (m_d_compute_get_tag(px) == m_exterior)
  {
    if (c->filter && c->filter->reject && c->filter->reject(c->filter, p))
    {
      p = m_d_compute_get_q(px);
      zp = m_d_compute_get_zq(px);
    }
  }
  int q1 = creal(zp) > 0;
  int q2 = (cimag(zp) > 0) != q1;
  int w = 0;
  int b = 0;
  if (m_d_compute_get_tag(px) == m_exterior)
  {
    double _Complex z = m_d_compute_get_z(px);
    w = cimag(z) > 0;
    int n = m_d_compute_get_n(px);
    b = n & 1;
  }
  return m_pixel_blacken
    ( m_pixel_blacken(m_pixel_whiten(m_pixel_hsva
        ( c->hue * (p - 1)
        , c->sat * (0.8 + 0.2 * q1)
        , c->val * (0.8 + 0.2 * q2)
        , 1
        ), 0.2 * w), 0.025 * b)
    , 1 - tanh(m_d_compute_get_de(px) / cabs(dc))
    );
}

extern m_d_colour_t *m_d_colour_domain(double hue, double sat, double val, m_period_filter_t *filter)
{
  m_d_colour_domain_t *c = malloc(sizeof(*c));
  c->colour.colour = m_d_colour_domain_do;
  c->hue = hue;
  c->sat = sat;
  c->val = val;
  c->filter = filter;
  return &c->colour;
}
