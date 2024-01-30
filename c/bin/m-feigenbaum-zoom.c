#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  const complex double c = -1.401155189093314712;
  const double d = 4.669201609102990671853203821578;
  double r = 8;
  int w = 512;
  int h = 512;
  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 600;
  int maxiters = 100000;
  m_image *image = m_image_new(w, h);
  if (image) {
    m_d_colour_t *colour = m_d_colour_minimal(red, black, white);
    if (colour) {
      for (int frame = 0; frame < 15; ++frame) {
        m_d_transform *rect = m_d_transform_rectangular(w, h, c, r);
        double o = -3.0 / 8.0;
        m_d_transform *move = m_d_transform_linear((w + I * h) * o, 1);
        m_d_transform *transform = m_d_transform_compose(move, rect);
        m_d_render_scanline(image, transform, er, maxiters, colour);
        char filename[100];
        snprintf(filename, 100, "%02d.png", frame);
        m_image_save_png(image, filename);
        m_d_transform_delete(transform);
        m_d_transform_delete(rect);
        m_d_transform_delete(move);
        r /= d;
      }
      m_d_colour_delete(colour);
    }
    m_image_delete(image);
  }
  return 0;
}
