/*
gcc -std=c99 -Wall -pedantic -Wextra -O3 -march=native -fopenmp \
    -o m-misiurewicz-domains m-misiurewicz-domains.c \
    `PKG_CONFIG_PATH=/path/to/mandelbrot/opt/lib/pkgconfig \
    pkg-config --cflags --libs mandelbrot-graphics`
*/
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mandelbrot-numerics.h>
#include <mandelbrot-graphics.h>

static inline double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

extern int main(int argc, char **argv) {
  if (argc != 9) {
    fprintf(stderr, "usage: %s out.png width height creal cimag radius maxiters period\n", argv[0]);
    return 1;
  }
  int width = atoi(argv[2]);
  int height = atoi(argv[3]);
  double _Complex center = atof(argv[4]) + I * atof(argv[5]);
  double radius = atof(argv[6]);
  int maxiters = atoi(argv[7]);
  int period = atoi(argv[8]);
  const char *filename = argv[1];

  int *pps = calloc(1, maxiters * sizeof(*pps));
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
      complex double zp = c;
      complex double dcp = 1;
      complex double dz = 0;
      int pp = 0;
      double mp2 = 1.0 / 0.0;
      double mz2 = 1.0 / 0.0;
      int np = 0;
      for (int n = 0; n < period; ++n) {
        double z2 = cabs2(z);
        if (z2 < mz2) {
          mz2 = z2;
          int p = n + 1;
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
      if (de == 0) {
        for (int n = 0; n < maxiters - period; ++n) {
          double z2 = cabs2(z);
          if (z2 < mz2) {
            mz2 = z2;
            int p = n + 1 + period;
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
          double p2 = cabs2(z - zp);
          if (p2 < mp2) {
            mp2 = p2;
            pp = n;
          }
          dc = 2 * z * dc + 1;
          z = z * z + c;
          dcp = 2 * zp * dcp + 1;
          zp = zp * zp + c;
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
      double hue = 0;
      if (pp > 0)
      {
        m_d_misiurewicz(&c, c, pp, period, 16);
        m_d_misiurewicz_naive(&c, c, pp, period, 16);
        z = 0;
        for (int n = 0; n < pp; ++n) z = z * z + c;
        complex double dz = 1;
        for (int n = 0; n < period; ++n)
        {
          dz = 2 * z * dz;
          z = z * z + c;
        }
        hue = fabs(carg(dz)) / log(cabs(dz));
      }
      double sat = bias == m_interior ? 0 : (pp > 0) * (0.25 + 0.75 /
        pow(1.0 + sqrt(mp2 / pixelspacing2) / 1024, 2.0));
      m_image_plot(img, i, j,
        m_pixel_hsva(hue, sat, tanh(de + 0.5), 1));
      #pragma omp atomic
      pps[pp]++;
    }
    free(zs);
    free(ps);
  }
  m_image_dirty(img);
  m_image_save_png(img, filename);
  m_d_transform_delete(transform);
  m_image_delete(img);
  for (int pp = 0; pp < maxiters; ++pp)
    if (pps[pp])
      printf("%d\t%d\n", pps[pp], pp);
  free(pps);
  return 0;
}
