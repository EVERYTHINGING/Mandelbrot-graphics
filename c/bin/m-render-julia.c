#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv) {
  if (argc != 10) {
    fprintf(stderr, "usage: %s out.png width height zreal zimag radius maxiters creal cimag\n", argv[0]);
    return 1;
  }
  const char *filename = argv[1];
  int w = atoi(argv[2]);
  int h = atoi(argv[3]);
  double _Complex z0 = atof(argv[4]) + I * atof(argv[5]);
  double r = atof(argv[6]);
  int maxiters = atoi(argv[7]);
  double _Complex c = atof(argv[8]) + I * atof(argv[9]);
  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 25;
  double er2 = er * er;
  int retval = 1;
  m_image *image = m_image_new(w, h);
  if (image) {
    m_d_transform *transform = m_d_transform_rectangular(w, h, z0, r);
    if (transform) {
      m_image_flush(image);
      #pragma omp parallel for
      for (int y = 0; y < h; ++y)
      {
        for (int x = 0; x < w; ++x)
        {
          double _Complex z = x + I * y;
          double _Complex dz = 1;
          m_d_transform_forward(transform, &z, &dz);
          for (int i = 0; i < maxiters; ++i)
          {
            if (cnorm(z) > er2)
              break;
            dz = 2 * z * dz;
            z = z * z + c;
          }
          m_pixel_t px = red;
          if (cnorm(z) > er2)
          {
            double de = tanh(cabs(z) * log(cabs(z)) / cabs(dz));
            px = m_pixel_mix(black, white, de);
          }
          m_image_plot(image, x, y, px);
        }
      }
      m_image_dirty(image);
      m_image_save_png(image, filename);
      m_d_transform_delete(transform);
      retval = 0;
    }
    m_image_delete(image);
  }
  return retval;
}
