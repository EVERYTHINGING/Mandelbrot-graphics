#include <mandelbrot-numerics.h>
#include <mandelbrot-graphics.h>
#include <stdio.h>

extern int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  int w = 1500;
  int h = 500;
  double r = 2.25;
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 100;
  int maxiters = 100100;
  const char *filename = "out.png";

  m_d_transform *rect = m_d_transform_rectangular(w, h, -2.5 * r, r);
  m_d_transform *exponential = m_d_transform_exponential(-1.401155189093314712);
  m_d_transform *transform = m_d_transform_compose(rect, exponential);
  m_d_colour_t *colour = m_d_colour_minimal(white, black, white);
  m_image *image = m_image_new(w, h);
  m_d_render_scanline(image, transform, er, maxiters, colour);
  m_image_save_png(image, filename);
  return 0;
}
