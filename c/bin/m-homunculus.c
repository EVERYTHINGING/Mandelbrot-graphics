#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;
  double _Complex c0 = -1.9409856638151786271684397e+00 + 6.4820395780451436662598436e-04 * I;
  double r0 = 5e-9;
  int w = 1920 * 4;
  int h = 1080 * 4;
  double er = 600;
  double maxiters = 100000;
  m_image *image = m_image_new(w, h);
  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  m_d_transform *transform = m_d_transform_rectangular(w, h, c0, r0);
  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    for (int i = 0; i < w; ++i)
    {
      double _Complex z = 0;
      double _Complex c = i + I * j;
      double _Complex dc0 = 1;
      m_d_transform_forward(transform, &c, &dc0);
      int p = 0;
      double _Complex z0 = 0;
      double mz1 = 1.0 / 0.0;
      double _Complex z1 = 0;
      bool e = false;
      double de = -1;
      for (int n = 1; n < maxiters; ++n)
      {
        z = z * z + c;
        double mz = cnorm(z);
        if (mz < mz1)
        {
          z0 = z1;
          mz1 = mz;
          z1 = z;
          p = n;
        }
        if (mz > er)
        {
          e = true;
          break;
        }
      }
      if (e)
      {
      double _Complex a = z1 / z0;
      double _Complex mu = c, tmp, b;
      m_d_nucleus(&mu, mu, p, 64);
      m_d_interior(&tmp, &b, mu, mu, -1, p, 64);
      b = 0.5 * (mu + b);
      double _Complex s = m_d_size(mu, p);
      double t = m_d_domain_size(mu, p);
      double f = 2 * cabs(s) / cabs(t);
      double g = pow(fmin(fmax(cabs(a), 0), 1), 64);
      double _Complex ch = (f * (c - b) + b) * (1 - g) + g * c;
      double _Complex dch = (f * (1 - g) + g) * dc0;
      z = 0;
      double _Complex dc = 0;
      for (int n = 1; n < maxiters; ++n)
      {
        dc = 2 * dc * z + 1;
        z = z * z + ch;
        double mz = cnorm(z);
        if (mz > er)
        {
          de = 2 * cabs(z) * log(cabs(z)) / cabs(dc * dch);
          break;
        }
      }
      }
      m_pixel_t colour = red;
      if (de >= 0)
      {
        colour = m_pixel_mix(black, white, tanh(de));
      }
      m_image_plot(image, i, j, colour);
    }
  }
  m_image_save_png(image, "m-homunculus.png");
  m_d_transform_delete(transform);
  m_image_delete(image);
  return 0;
}
