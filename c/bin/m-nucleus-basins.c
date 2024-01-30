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
      m_d_nucleus(&z, c, period, 64);
      double hue = carg(z - center) / 6.283185307179586; hue -= floor(hue);
      double sat = cabs(z - c) < radius / 16 ? 1 : 0.25;
      double val = tanh(de + 1);
      m_image_plot(img, i, j, m_pixel_hsva(hue, sat, val, 1));
    }
    free(zs);
    free(ps);
  }
  m_image_dirty(img);
  m_image_save_png(img, filename);
  m_d_transform_delete(transform);
  m_image_delete(img);
  return 0;
}
