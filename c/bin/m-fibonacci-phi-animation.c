#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  const double twopi = 6.283185307179586;
  const double phi = 1.618033988749895;
  complex double u = cexp(I * twopi / phi) / 2;
  complex double c = u * (1 - u);
  double f = 1 / (phi * phi);
  f *= f;
  int w = 512;
  int h = 512;
  double r = 4;
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 600;
  int maxiters = 1000000;
  m_image *image = m_image_new(w, h);
  if (image) {
    m_d_colour_t *colour = m_d_colour_minimal(white, black, white);
    if (colour) {
      for (int i = 0; i < 8; ++i) {
        m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);
        if (transform) {
          m_d_render_scanline(image, transform, er, maxiters, colour);
          char filename[100];
          snprintf(filename, 100, "%03d.png", i);
          m_image_save_png(image, filename);
          m_d_transform_delete(transform);
        }
        r *= f;
      }
      m_d_colour_delete(colour);
    }
    m_image_delete(image);
  }
  return 0;
}
