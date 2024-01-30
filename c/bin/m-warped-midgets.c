#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <mandelbrot-symbolics.h>
#include <mandelbrot-numerics.h>
#include <mandelbrot-graphics.h>

int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;
  int w = 256;
  int h = 256;
  int sharpness = 4;
  double er = 600;
  int maxiters = 1 << 18;
  m_image *image = m_image_new(w, h);
  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  mpq_t angle;
  mpq_init(angle);
  m_binangle bangle;
  m_binangle_init(&bangle);
  if (image)
  {
    cairo_surface_t *surface = m_image_surface(image);
    cairo_t *cr = cairo_create(surface);
    cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 24);
    cairo_set_source_rgba(cr, 1, 0, 0, 1);
    m_d_colour_t *colour = m_d_colour_minimal(red, black, white);
    if (colour) {
      for (int period = 4; period < 1 << 14; period <<= 1)
      {
        mpz_set_ui(bangle.pre.bits, 0); bangle.pre.length = 0;
        mpz_set_ui(bangle.per.bits, 3); bangle.per.length = period;
        m_binangle_to_rational(angle, &bangle);
        complex double ray = m_d_exray_in_do(angle, sharpness, 8 * sharpness * period, 64);
        complex double mc = 0;
        m_d_nucleus(&mc, ray, period, 64);
        double r = cabs(m_d_size(mc, period)) * 2;
        printf("%7d %.18e + %.18e i @ %.3e\n", period, creal(mc), cimag(mc), r);
        m_d_transform *transform = m_d_transform_rectangular(w, h, mc, r);
        if (transform)
        {
          m_d_render_scanline(image, transform, er, maxiters, colour);
          m_image_dirty(image);
          cairo_move_to(cr, 24, 24);
          char text[100];
          snprintf(text, 100, "%d", period);
          cairo_show_text(cr, text);
          cairo_fill(cr);
          m_image_flush(image);
          char filename[100];
          snprintf(filename, 100, "%06d.png", period);
          m_image_save_png(image, filename);
          m_d_transform_delete(transform);
        }
      }
      m_d_colour_delete(colour);
    }
    m_image_delete(image);
  }
  m_binangle_clear(&bangle);
  mpq_clear(angle);
  return 0;
}
