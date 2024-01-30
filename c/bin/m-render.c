#include <stdio.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv) {
  if (argc != 9) {
    fprintf(stderr, "usage: %s out.png width height creal cimag radius maxiters rgb\n", argv[0]);
    return 1;
  }
  const char *filename = argv[1];
  int w = atoi(argv[2]);
  int h = atoi(argv[3]);
  double _Complex c = atof(argv[4]) + I * atof(argv[5]);
  double r = atof(argv[6]);
  int maxiters = atoi(argv[7]);
  int rgb = atoi(argv[8]);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double phi = (sqrt(5) + 1) / 2;
  double gold = 1 / (phi * phi);
  double er = 25;
  int retval = 1;
  m_image *image = m_image_new(w, h);
  if (image) {
    m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);
    if (transform) {
      m_d_colour_t *colour;
      if (rgb)
        colour = m_d_colour_domain(gold / rgb, 0.5, 1, 0);
      else
        colour = m_d_colour_minimal(white, black, white);
      if (colour) {
        m_d_render_scanline(image, transform, er, maxiters, colour);
        m_image_save_png(image, filename);
        retval = 0;
        m_d_colour_delete(colour);
      }
      m_d_transform_delete(transform);
    }
    m_image_delete(image);
  }
  return retval;
}
