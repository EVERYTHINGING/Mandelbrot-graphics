#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  complex double c = -0.75 + I * 0.0;
  double r = 1.5;
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 600;
  int maxiters = 1000000;
  for (int b = 8; b < 16; ++b)
  {
    int w = 1 << b;
    int h = 1 << b;
    m_image *image = m_image_new(w, h);
    if (image) {
      m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);
      if (transform) {
        complex double c0 = 0;
        complex double dc0 = 1;
        m_d_transform_forward(transform, &c0, &dc0);
        double px = cabs(dc0);
        m_d_colour_t *colour = m_d_colour_minimal(black, black, white);
        if (colour) {
          m_d_render_scanline(image, transform, er, maxiters, colour);
          char filename[100];
          snprintf(filename, 100, "%d.png", b);
          m_image_save_png(image, filename);
          double s = 0;
          double n = 0;
          #pragma omp parallel for reduction(+:s) reduction(+:n) schedule(static, 1)
          for (int y1 = 0; y1 < h; ++y1)
          {
            double s0 = 0;
            double s1 = 0;
            for (int x1 = 0; x1 < w; ++x1)
              if (white != m_image_peek(image, x1, y1))
                for (int y2 = 0; y2 < h; ++y2)
                  for (int x2 = 0; x2 < w; ++x2)
                    if (white != m_image_peek(image, x2, y2))
                    {
                      s0 += 1.0;
                      s1 += hypot(x1 - x2, y1 - y2);
                    }
            s1 *= px;
            n = n + s0;
            s = s + s1;
          }
          s /= n;
          printf("%d\t%.18f\n", b, s);
          fflush(stdout);
          m_d_colour_delete(colour);
        }
        m_d_transform_delete(transform);
      }
      m_image_delete(image);
    }
  }
  return 0;
}
