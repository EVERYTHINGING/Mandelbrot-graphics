#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;
  int w = 256;
  int h = 256;
  double er = 600;
  double maxiters = 100000;
  m_image *image = m_image_new(w, h);
  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  m_d_colour_t *colour = m_d_colour_minimal(red, black, white);
  m_d_transform *rect = m_d_transform_rectangular(w, h, 0, 2);
  m_d_transform *card = m_d_transform_cardioid();
  m_d_transform *uncard = m_d_transform_cardioid();
  m_d_transform_invert(uncard);
  m_d_transform *line = m_d_transform_moebius3(1, -I, -1);
  m_d_transform *unline = m_d_transform_moebius3(1, -I, -1);
  m_d_transform_invert(unline);
  for (int x = 0; x < 256; ++x)
  {
    m_d_transform *translate = m_d_transform_linear(x / 64.0, 1.0);
    m_d_transform *t0 = m_d_transform_compose(rect, uncard);
    m_d_transform *t1 = m_d_transform_compose(t0, unline);
    m_d_transform *t2 = m_d_transform_compose(t1, translate);
    m_d_transform *t3 = m_d_transform_compose(t2, line);
    m_d_transform *t4 = m_d_transform_compose(t3, card);
    m_d_render_scanline(image, t4, er, maxiters, colour);
    char filename[100];
    snprintf(filename, 100, "%03d.png", x);
    m_image_save_png(image, filename);
    m_d_transform_delete(translate);
    m_d_transform_delete(t0);
    m_d_transform_delete(t1);
    m_d_transform_delete(t2);
    m_d_transform_delete(t3);
    m_d_transform_delete(t4);
  }
  m_image_delete(image);
  m_d_colour_delete(colour);
  m_d_transform_delete(rect);
  m_d_transform_delete(line);
  m_d_transform_delete(unline);
  m_d_transform_delete(card);
  m_d_transform_delete(uncard);
  return 0;
}
