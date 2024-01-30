#include <mandelbrot-graphics.h>
#include <mandelbrot-numerics.h>
#include "m_d_util.h"

struct m_d_partial {
  double _Complex z;
  int p;
};
typedef struct m_d_partial m_d_partial;

struct m_d_compute {
  m_compute_t tag, bias;
  m_period_filter_t *filter;
  m_d_partial *partials;
  int npartials, np, n, p, q;
  double er2, mz2, mzq2, de;
  double _Complex c, z, dc, dz, zp, zq;
};

extern m_d_compute *m_d_compute_alloc(int npartials) {
  m_d_compute *px = malloc(sizeof(*px));
  if (! px) {
    return 0;
  }
  if (npartials > 0) {
    px->partials = malloc(npartials * sizeof(*(px->partials)));
    px->npartials = npartials;
    if (! px->partials) {
      free(px);
      return 0;
    }
  } else {
    px->partials = 0;
    px->npartials = 0;
  }
  px->tag = m_unknown;
  return px;
}

extern void m_d_compute_free(m_d_compute *px) {
  if (px) {
    if (px->partials) {
      free(px->partials);
    }
    free(px);
  }
}

extern void m_d_compute_init(m_d_compute *px, m_compute_t bias, double er, double _Complex c, m_period_filter_t *filter) {
  if (! px) { return; }
  px->tag = m_unknown;
  px->bias = bias;
  px->er2 = er * er;
  px->filter = filter;
  px->mz2 = 1.0 / 0.0;
  px->mzq2 = 1.0 / 0.0;
  px->c = c;
  px->z = 0;
  px->dc = 0;
  px->dz = 0;
  px->zp = 0;
  px->zq = 0;
  px->n = 0;
  px->p = 0;
  px->q = 0;
  px->np = 0;
  px->de = -1;
}

extern void m_d_compute_clear(m_d_compute *px) {
  px->tag = m_unknown;
  px->np = 0;
}

extern bool m_d_compute_step(m_d_compute *px, int steps) {
  if (! px) {
    return false;
  }
  if (px->tag != m_unknown) {
    return true;
  }
  double er2 = px->er2;
  double _Complex c = px->c;
  double _Complex z = px->z;
  double _Complex dc = px->dc;
  double _Complex zp = px->zp;
  double _Complex zq = px->zq;
  double mz2 = px->mz2;
  double mzq2 = px->mzq2;
  int p = px->p;
  int q = px->q;
  for (int i = 1; i <= steps; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c;
    double z2 = cabs2(z);
    if (z2 < mzq2 && px->filter && px->filter->accept && px->filter->accept(px->filter, px->n + i))
    {
      mzq2 = z2;
      q = px->n + i;
      zq = z;
    }
    if (z2 < mz2) {
      double atom_domain_radius_squared = z2 / mz2;
      mz2 = z2;
      p = px->n + i;
      zp = z;
      if (atom_domain_radius_squared <= 0.25) {
        if (px->bias == m_interior) {
          double _Complex dz = 0;
          double de = -1;
          if (m_d_interior_de(&de, &dz, z, c, p, 64)) {
            px->tag = m_interior;
            px->p = p;
            px->z = z;
            px->dz = dz;
            px->zp = zp;
            px->de = de;
            return true;
          }
        } else {
          if (px->partials && px->np < px->npartials) {
            px->partials[px->np].z = z;
            px->partials[px->np].p = p;
            px->np = px->np + 1;
          }
        }
      }
    }
    if (! (z2 < er2)) {
      px->tag = m_exterior;
      px->n = px->n + i;
      px->p = p;
      px->q = q;
      px->z = z;
      px->zp = zp;
      px->zq = zq;
      px->dc = dc;
      px->de = 2 * cabs(z) * log(cabs(z)) / cabs(dc);
      return true;
    }
  }
  if (px->bias != m_interior && px->partials) {
    for (int i = 0; i < px->np; ++i) {
      z = px->partials[i].z;
      zp = z;
      int p = px->partials[i].p;
      double _Complex dz = 0;
      double de = -1;
      if (m_d_interior_de(&de, &dz, z, c, p, 64)) {
        px->tag = m_interior;
        px->p = p;
        px->z = z;
        px->dz = dz;
        px->zp = zp;
        px->de = de;
        return true;
      }
    }
  }
  px->tag = m_unknown;
  px->n = px->n + steps;
  px->p = p;
  px->q = q;
  px->mz2 = mz2;
  px->mzq2 = mzq2;
  px->z = z;
  px->dc = dc;
  px->zp = zp;
  px->zq = zq;
  return false;
}

extern m_compute_t m_d_compute_get_tag(const m_d_compute *px) {
  if (! px) { return m_unknown; }
  return px->tag;
}

extern int m_d_compute_get_n(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->n;
}

extern int m_d_compute_get_p(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->p;
}

extern int m_d_compute_get_q(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->q;
}

extern double _Complex m_d_compute_get_c(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->c;
}

extern double _Complex m_d_compute_get_z(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->z;
}

extern double _Complex m_d_compute_get_dc(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->dc;
}

extern double _Complex m_d_compute_get_dz(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->dz;
}

extern double _Complex m_d_compute_get_zp(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->zp;
}

extern double _Complex m_d_compute_get_zq(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->zq;
}

extern double m_d_compute_get_de(const m_d_compute *px) {
  if (! px) { return 0; }
  return px->de;
}
