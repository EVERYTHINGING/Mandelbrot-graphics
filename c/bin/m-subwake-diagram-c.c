#include <mandelbrot-graphics.h>
#include <mandelbrot-numerics.h>
#include <mandelbrot-symbolics.h>
#include <cairo.h>
#include <stdio.h>

const double twopi = 6.283185307179586;

void draw_label(m_image *image, m_d_transform *transform, double _Complex c0, const char *text, double pt, m_pixel_t colour, double angle) {
  double _Complex c = c0;
  double _Complex dc = 1;
  m_d_transform_reverse(transform, &c, &dc);
  cairo_surface_t *surface = m_image_surface(image);
  cairo_t *cr = cairo_create(surface);
  cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, pt);
  cairo_text_extents_t te;
  cairo_text_extents(cr, text, &te);
  cairo_translate(cr, creal(c), cimag(c));
  cairo_rotate(cr, angle);
  cairo_translate(cr, - te.x_bearing - te.width / 2, - te.y_bearing - te.height / 2);
  cairo_text_path(cr, text);
  cairo_set_source_rgba(cr, m_pixel_red(colour), m_pixel_green(colour), m_pixel_blue(colour), m_pixel_alpha(colour));
  cairo_fill(cr);
  cairo_destroy(cr);
}

void draw_external_ray(m_image *image, m_d_transform *transform, const char *angle, m_pixel_t colour, double pt, m_pixel_t lcolour, bool plabel, m_pixel_t acolour, const char *addr) {
  int maxiters = 1024;

  m_binangle btheta;
  m_binangle_init(&btheta);
  m_binangle_from_string(&btheta, angle);
  int period = btheta.per.length;

  mpq_t qtheta;
  mpq_init(qtheta);
  m_binangle_to_rational(qtheta, &btheta);
  m_binangle_clear(&btheta);
  m_d_exray_in *ray = m_d_exray_in_new(qtheta, 8);
  mpq_clear(qtheta);

  cairo_surface_t *surface = m_image_surface(image);
  cairo_t *cr = cairo_create(surface);
  cairo_set_source_rgba(cr, m_pixel_red(colour), m_pixel_green(colour), m_pixel_blue(colour), m_pixel_alpha(colour));
  bool first = true;
  double _Complex c1 = 0;
  for (int i = 0; i < maxiters; ++i) {
    if (m_failed == m_d_exray_in_step(ray, 64)) {
      break;
    }
    double _Complex c = m_d_exray_in_get(ray);
    c1 = c;
    double _Complex dc = 1;
    double t = cimag(c) > 0 ? twopi / 4 : -twopi / 4;
    m_d_transform_reverse(transform, &c, &dc);
    double h = fabs(cimag(c) - m_image_get_height(image) / 2);
    if (h > m_image_get_height(image) / 6) {
      continue;
    }
    if (first) {
      cairo_save(cr);
      cairo_translate(cr, creal(c), cimag(c));
      cairo_rotate(cr, -t);
      cairo_select_font_face(cr, "LMMono10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
      cairo_set_font_size(cr, 48);
      cairo_text_path(cr, angle);
      cairo_fill(cr);
      cairo_restore(cr);
      cairo_move_to(cr, creal(c), cimag(c));
      first = false;
    } else {
      cairo_line_to(cr, creal(c), cimag(c));
    }
  }
  cairo_stroke(cr);
  cairo_destroy(cr);

  if (plabel) {
    double _Complex nucleus;
    m_d_nucleus(&nucleus, c1, period, 64);  
    char speriod[100];
    snprintf(speriod, 100, "%d", period);
    draw_label(image, transform, nucleus, speriod, pt, lcolour, 0);
    draw_label(image, transform, nucleus + 0.005 + I * 0.005, addr, pt / 2, acolour, -twopi / 12);
  }
}

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  int w = 16384;
  int h = 1024;

  double _Complex c3, c4, c5, c6, c7, c8;
  m_d_nucleus(&c3, -2, 3, 64);
  m_d_nucleus(&c4, -2, 4, 64);
  m_d_nucleus(&c5, -2, 5, 64);
  m_d_nucleus(&c6, -2, 6, 64);
  m_d_nucleus(&c7, -2, 7, 64);
  m_d_nucleus(&c8, -1.4, 8, 64);

  double er = 600;
  int maxiters = 8192;
  const char *filename = "subwake-diagram-c.png";

  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t green = m_pixel_rgba(0, 0.5, 0, 1);
  m_pixel_t blue  = m_pixel_rgba(0, 0, 1, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double pt = 96;

  int retval = 1;
  m_image *image = m_image_new(w, h);
  if (image) {
    double _Complex c = (c7 + c8) / 2;
    double r = (2.03 + creal(c)) * h / (double) w;
    m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);

    if (transform) {
      m_d_colour_t *colour = m_d_colour_minimal(white, black, white);
      if (colour) {
        m_d_render_scanline(image, transform, er, maxiters, colour);
        /*
        $ cat feature-database.csv | grep True | grep ^[1-7], | grep -v / | (
          IFS=,
          while read p i a n1 d1 n2 d2 rest
          do
            echo ", { \"`m-binangle-from-rational $n1/$d1`\", \"`m-binangle-from-rational $n2/$d2`\", \"$a\" }"
          done
        )
        */
        const char *angles[19][3] =
{ { ".(011)", ".(100)", "1 2 3" }
, { ".(0111)", ".(1000)", "1 2 3 4" }
, { ".(01111)", ".(10000)", "1 2 3 4 5" }
, { ".(01110)", ".(10001)", "1 2 3 5" }
, { ".(01101)", ".(10010)", "1 2 4 5" }
, { ".(011111)", ".(100000)", "1 2 3 4 5 6" }
, { ".(011110)", ".(100001)", "1 2 3 4 6" }
, { ".(011101)", ".(100010)", "1 2 3 5 6" }
, { ".(011010)", ".(100101)", "1 2 4 6" }
, { ".(0111111)", ".(1000000)", "1 2 3 4 5 6 7" }
, { ".(0111110)", ".(1000001)", "1 2 3 4 5 7" }
, { ".(0111101)", ".(1000010)", "1 2 3 4 6 7" }
, { ".(0111100)", ".(1000011)", "1 2 3 4 7" }
, { ".(0111011)", ".(1000100)", "1 2 3 5 6 7" }
, { ".(0111010)", ".(1000101)", "1 2 3 5 7" }
, { ".(0111001)", ".(1000110)", "1 2 3 6 7" }
, { ".(0110110)", ".(1001001)", "1 2 4 5 7" }
, { ".(0110101)", ".(1001010)", "1 2 4 6 7" }
, { ".(01101001)", ".(10010110)", "1 2 4 8" }
};
        for (int a = 0; a < 19; ++a) {
          for (int b = 0; b < 2; ++b) {
            draw_external_ray(image, transform, angles[a][b], red, pt, blue, b == 0, green, angles[a][2]);
          }
        }

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
