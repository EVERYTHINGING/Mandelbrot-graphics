#include <mandelbrot-graphics.h>
#include <mandelbrot-numerics.h>
#include <mandelbrot-symbolics.h>
#include <cairo.h>

const double twopi = 6.283185307179586;

void draw_label(m_image *image, m_d_transform *transform, double _Complex c0, const char *text, double pt, m_pixel_t colour) {
  double _Complex c = c0;
  double _Complex dc = 1;
  m_d_transform_reverse(transform, &c, &dc);
  cairo_surface_t *surface = m_image_surface(image);
  cairo_t *cr = cairo_create(surface);
  cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, pt);
  cairo_text_extents_t te;
  cairo_text_extents(cr, text, &te);
  cairo_move_to(cr, creal(c) - te.x_bearing - te.width / 2, cimag(c) - te.y_bearing - te.height / 2);
  cairo_text_path(cr, text);
  cairo_set_source_rgba(cr, m_pixel_red(colour), m_pixel_green(colour), m_pixel_blue(colour), m_pixel_alpha(colour));
  cairo_fill(cr);
  cairo_destroy(cr);
}

void draw_internal_ray(m_image *image, m_d_transform *transform, int period, double _Complex nucleus, const char *angle, double pt, m_pixel_t colour) {
  int steps = 128;
  mpq_t theta;
  mpq_init(theta);
  mpq_set_str(theta, angle, 10);
  mpq_canonicalize(theta);
  double a = twopi * mpq_get_d(theta);
  mpq_clear(theta);
  double _Complex interior = cos(a) + I * sin(a);

  double _Complex cl = 0, cl2 = 0;
  double _Complex c = nucleus;
  double _Complex z = c;
  cairo_surface_t *surface = m_image_surface(image);
  cairo_t *cr = cairo_create(surface);
  cairo_set_source_rgba(cr, m_pixel_red(colour), m_pixel_green(colour), m_pixel_blue(colour), m_pixel_alpha(colour));
  for (int i = 0; i < steps; ++i) {
    if (2 * i == steps) {
      cl = c;
    }
    if (2 * i == steps + 2) {
      cl2 = c;
    }
    double radius = (i + 0.5) / steps;
    m_d_interior(&z, &c, z, c, radius * interior, period, 64);
    double _Complex pc = c;
    double _Complex pdc = 1;
    m_d_transform_reverse(transform, &pc, &pdc);
    if (i == 0) {
      cairo_move_to(cr, creal(pc), cimag(pc));
    } else {
      cairo_line_to(cr, creal(pc), cimag(pc));
    }
  }
  cairo_stroke(cr);
  if (a != 0) {
    double t = carg(cl2 - cl);
    cairo_save(cr);
    double _Complex dcl = 1;
    m_d_transform_reverse(transform, &cl, &dcl);
    cairo_translate(cr, creal(cl), cimag(cl));
    cairo_rotate(cr, -t);
    cairo_translate(cr, 0, -pt/3);
    cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, pt);
    cairo_text_path(cr, angle);
    cairo_fill(cr);
    cairo_restore(cr);
  }
  cairo_destroy(cr);
}

void draw_external_ray(m_image *image, m_d_transform *transform, const char *angle, m_pixel_t colour, double dx, double dy) {
  int maxiters = 1024;
  double r = sqrt(2);

  m_binangle btheta;
  m_binangle_init(&btheta);
  m_binangle_from_string(&btheta, angle);
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
  for (int i = 0; i < maxiters; ++i) {
    if (m_failed == m_d_exray_in_step(ray, 64)) {
      break;
    }
    double _Complex c = m_d_exray_in_get(ray);
    if (cabs(c + 0.75) > r) {
      continue;
    }
    double t = carg(c + 0.75);
    double _Complex dc = 1;
    m_d_transform_reverse(transform, &c, &dc);
    if (first) {
      cairo_save(cr);
      cairo_translate(cr, creal(c) + dx, cimag(c) + dy);
      cairo_rotate(cr, -t);
      cairo_select_font_face(cr, "LMMono10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
      cairo_set_font_size(cr, 48);
      cairo_text_path(cr, angle);
      cairo_fill(cr);
      cairo_restore(cr);
      cairo_move_to(cr, creal(c) + dx, cimag(c) + dy);
      first = false;
    } else {
      cairo_line_to(cr, creal(c), cimag(c));
    }
  }
  cairo_stroke(cr);
  cairo_destroy(cr);
}

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  int w = 4096;
  int h = 4096;
  complex double c = -0.75;
  double r = 1.75;
  double er = 600;
  int maxiters = 8192;
  const char *filename = "subwake-diagram-a.png";

  m_pixel_t red   = m_pixel_rgba(1, 0, 0, 1);
  m_pixel_t green = m_pixel_rgba(0, 0.5, 0, 1);
  m_pixel_t blue  = m_pixel_rgba(0, 0, 1, 1);
  m_pixel_t black = m_pixel_rgba(0, 0, 0, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);

  int retval = 1;
  m_image *image = m_image_new(w, h);
  if (image) {
    m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);
    if (transform) {
      m_d_colour_t *colour = m_d_colour_minimal(white, black, white);
      if (colour) {
        m_d_render_scanline(image, transform, er, maxiters, colour);

        double _Complex c3, c4a, c4b, c5, c3c2, c2c3;
        m_d_nucleus(&c3, 0 + I * 1, 3, 64);
        m_d_nucleus(&c4a, 0.25 + 0.5 * I, 4, 64);
        m_d_nucleus(&c4b, 0.25 - 0.5 * I, 4, 64);
        m_d_nucleus(&c5,  0.3 + 0.3 * I, 5, 64);
        m_d_nucleus(&c3c2, c3 + I * 0.1, 6, 64);
        m_d_nucleus(&c2c3, -1 - 0.25 + 0.25 * I, 6, 64);

        double pt = 48;
        
        draw_internal_ray(image, transform, 1, 0, "1/2", pt, green);
        draw_internal_ray(image, transform, 1, 0, "1/3", pt, green);
        draw_internal_ray(image, transform, 1, 0, "1/4", pt, green);
        draw_internal_ray(image, transform, 1, 0, "1/5", pt, green);
        draw_internal_ray(image, transform, 1, 0, "3/4", pt, green);
        draw_internal_ray(image, transform, 2, -1, "0/1", pt, green);
        draw_internal_ray(image, transform, 2, -1, "1/3", pt, green);
        draw_internal_ray(image, transform, 3, c3, "0/1", 0.7 * pt, green);
        draw_internal_ray(image, transform, 3, c3, "1/2", 0.7 * pt, green);
        draw_internal_ray(image, transform, 3, c3, "1/3", 0.7 * pt, green);
        draw_internal_ray(image, transform, 3, c3, "1/4", 0.7 * pt, green);
        draw_internal_ray(image, transform, 3, c3, "3/4", 0.7 * pt, green);

        draw_external_ray(image, transform, ".(01)", red, 0, 0);
        draw_external_ray(image, transform, ".(10)", red, 0, 0);
        draw_external_ray(image, transform, ".(001)", red, 32, 32);
        draw_external_ray(image, transform, ".(010)", red, -48, 0);
        draw_external_ray(image, transform, ".(011)", red, 0, 0);
        draw_external_ray(image, transform, ".(100)", red, 0, 0);
        draw_external_ray(image, transform, ".(0001)", red, 0, -16);
        draw_external_ray(image, transform, ".(0010)", red, 16, 16);
        draw_external_ray(image, transform, ".(1101)", red, 0, 0);
        draw_external_ray(image, transform, ".(1110)", red, 0, 0);
        draw_external_ray(image, transform, ".(00001)", red, 0, 0);
        draw_external_ray(image, transform, ".(00010)", red, 0, 16);
        draw_external_ray(image, transform, ".(001010)", red, -32, -32);
        draw_external_ray(image, transform, ".(010001)", red, 0, 0);
        draw_external_ray(image, transform, ".(010110)", red, 0, 0);
        draw_external_ray(image, transform, ".(011001)", red, 0, 0);
        draw_external_ray(image, transform, ".(001001010)", red, -48, -48);
        draw_external_ray(image, transform, ".(001010001)", red, 0, 0);
        draw_external_ray(image, transform, ".(001001001010)", red, 0, 0);
        draw_external_ray(image, transform, ".(001001010001)", red, -16, -16);
        draw_external_ray(image, transform, ".(010010001010)", red, 32, 0);
        draw_external_ray(image, transform, ".(010010010001)", red, 0, 0);

        draw_label(image, transform,  0,   "1", 6 * pt, blue);
        draw_label(image, transform, -1,   "2", 3 * pt, blue);
        draw_label(image, transform, c3,   "3", 2 * pt, blue);
        draw_label(image, transform, c4a,  "4", 1.5 * pt, blue);
        draw_label(image, transform, c4b,  "4", 1.5 * pt, blue);
        draw_label(image, transform, c5,   "5", pt, blue);
        draw_label(image, transform, c2c3, "6", pt, blue);
        draw_label(image, transform, c3c2, "6", pt, blue);

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
