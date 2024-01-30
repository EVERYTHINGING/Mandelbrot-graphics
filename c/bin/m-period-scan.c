#include <stdio.h>
#include <mandelbrot-graphics.h>
#include <mandelbrot-numerics.h>

struct atom
{
  double _Complex nucleus;
  double size;
  double domain_size;
  int period;
  m_shape shape;
};

int cmp_atom_size(const void *a, const void *b)
{
  const struct atom *x = a;
  const struct atom *y = b;
  double s = x->size;
  double t = y->size;
  if (s > t) return -1;
  if (s < t) return 1;
  return 0;
}

int cmp_atom_domain_size(const void *a, const void *b)
{
  const struct atom *x = a;
  const struct atom *y = b;
  double s = x->domain_size;
  double t = y->domain_size;
  if (s > t) return -1;
  if (s < t) return 1;
  return 0;
}

int cmp_atom_period(const void *a, const void *b)
{
  const struct atom *x = a;
  const struct atom *y = b;
  int s = x->period;
  int t = y->period;
  if (s > t) return 1;
  if (s < t) return -1;
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 14) {
    fprintf(stderr, "usage: %s out.png width height creal cimag radius maxiters \\\n"
    "  mingridsize minfontsize maxfontsize maxatoms periodmod periodneq\n", argv[0]);
    return 1;
  }
  const char *filename = argv[1];
  int w = atoi(argv[2]);
  int h = atoi(argv[3]);
  double _Complex c = atof(argv[4]) + I * atof(argv[5]);
  double r = atof(argv[6]);
  int maxiters = atoi(argv[7]);
  int mingridsize = atoi(argv[8]);
  double minfontsize = atof(argv[9]);
  double maxfontsize = atof(argv[10]);
  int maxatoms = atoi(argv[11]);
  int periodmod = atoi(argv[12]);
  int periodneq = atoi(argv[13]);
  m_pixel_t grey = m_pixel_rgba(0.75, 0.75, 0.75, 1);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);
  double er = 600;
  m_image *image = m_image_new(w, h);
  if (image)
  {
    m_d_transform *transform = m_d_transform_rectangular(w, h, c, r);
    if (transform) {
      m_d_colour_t *colour = m_d_colour_minimal(white, grey, white);
      if (colour) {
        // render image
        m_d_render_scanline(image, transform, er, maxiters, colour);
        m_image_dirty(image);
        // scan for periods
        int atoms = 0;
        for (int grid = mingridsize << 8; grid >= mingridsize; grid >>= 1)
          for (int y = grid/2; y < h; y += grid)
            for (int x = grid/2; x < w; x += grid)
              atoms++;
        struct atom *as = calloc(1, atoms * sizeof(*as));
        atoms = 0;
        for (int grid = mingridsize << 8; grid >= mingridsize; grid >>= 1)
          for (int y = grid/2; y < h; y += grid)
            for (int x = grid/2; x < w; x += grid)
            {
              double _Complex c0 = x + I * y;
              double _Complex dc0 = grid;
              m_d_transform_forward(transform, &c0, &dc0);
              int p = m_d_box_period_do(c0, 4.0 * cabs(dc0), maxiters);
              if (p > 0)
                if (m_converged == m_d_nucleus(&c0, c0, p, 16))
                {
                  as[atoms].period = m_d_box_period_do(c0, 0.001 * cabs(dc0), 2 * p);
                  if (as[atoms].period > 0)
                  {
                    as[atoms].nucleus = c0;
                    as[atoms].size = cabs(m_d_size(c0, as[atoms].period));
                    as[atoms].domain_size = m_d_domain_size(c0, as[atoms].period);
                    as[atoms].shape = m_d_shape_discriminant(m_d_shape_estimate(c0, as[atoms].period));
                    atoms++;
                  }
                }
            }
        qsort(as, atoms, sizeof(*as), cmp_atom_domain_size);
        // prepare deduplication buffer
        m_image *maskimage = m_image_new(w, h);
        cairo_surface_t *masksurface = m_image_surface(maskimage);
        cairo_t *mask = cairo_create(masksurface);
        cairo_set_source_rgba(mask, 1, 1, 1, 1);
        cairo_paint(mask);
        cairo_set_source_rgba(mask, 0, 0, 0, 1);
        // prepare image
        cairo_surface_t *surface = m_image_surface(image);
        cairo_t *cr = cairo_create(surface);
        cairo_set_source_rgba(cr, 0.5, 0.25, 0.25, 0.75);
        cairo_set_operator(cr, CAIRO_OPERATOR_MULTIPLY);
        cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        int printed = 0;
        for (int a = 0; printed < maxatoms && a < atoms; ++a)
        {
          // convert to pixel coordinates
          int p = as[a].period;
          double _Complex c0 = as[a].nucleus;
          double _Complex dc0 = p == 1 ? 1 : as[a].domain_size;
          m_d_transform_reverse(transform, &c0, &dc0);
          double x = creal(c0);
          double y = cimag(c0);
          // calculate text size
          double fs = (as[a].shape == m_cardioid ? 1 : 0.5) * cabs(dc0);
          if (periodneq >= 0 && (p % periodmod) != periodneq)
            fs = 8 * log2(fs) + maxfontsize;
          fs = fmax(fs, minfontsize);
          cairo_set_font_size(cr, fs);
          // prevent drawing duplicates
          m_image_flush(maskimage);
          if (! (m_image_in_bounds(maskimage, x, y) && white == m_image_peek(maskimage, x, y)))
            continue;
          cairo_arc(mask, x, y, fs, 0, 6.283185307179586);
          cairo_close_path(mask);
          cairo_fill(mask);
          // draw text centered on point
          char sp[100];
          snprintf(sp, 100, "%d", p);
          cairo_text_extents_t e;
          cairo_text_extents(cr, sp, &e);
          x = x - e.width / 2.0;
          y = y - e.height / 2.0;
          cairo_move_to(cr, x - e.x_bearing, y - e.y_bearing);
          cairo_show_text(cr, sp);
          cairo_fill(cr);
          printed++;
        }
        // cleanup
        cairo_destroy(mask);
        m_image_delete(maskimage);
        cairo_destroy(cr);
        m_image_flush(image);
        m_image_save_png(image, filename);
        m_d_colour_delete(colour);
      }
      m_d_transform_delete(transform);
    }
    m_image_delete(image);
  }
  return 0;
}
