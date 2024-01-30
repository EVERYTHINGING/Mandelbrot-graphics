#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mandelbrot-graphics.h>
#include <mandelbrot-numerics.h>

#define IMAGE_E

void initialize_cs(int m, int n, m_d_transform *t, double _Complex *cs)
{
  #pragma omp parallel for
  for (int j = 0; j < n; ++j)
  {
    for (int i = 0; i < m; ++i)
    {
      double _Complex c = i + I * j;
      double _Complex dc = 1;
      m_d_transform_forward(t, &c, &dc);
      int k = i + j * m;
      cs[k] = c;
    }
  }
}

void step_zs(int mn, char *qs, double _Complex *zs, const double _Complex *cs)
{
  #pragma omp parallel for
  for (int i = 0; i < mn; ++i)
  {
    // load
    double _Complex c = cs[i];
    double _Complex z = zs[i];
    // step
    z = z * z + c;
    // compute quadrant
    char q = 1 << ((creal(z) > 0) | ((cimag(z) > 0) << 1));
    // store
    zs[i] = z;
    qs[i] = q;
  }
}

int scan_for_zeroes(int m, int n, int ip, int *ops, double _Complex *ocs, const char *qs, const double _Complex *zs, const double _Complex *ics)
{
  int o = 0;
  // loop over image interior, to avoid tests in inner 3x3 loop
  #pragma omp parallel for
  for (int j = 1; j < n - 1; ++j)
  {
    for (int i = 1; i < m - 1; ++i)
    {
      // find where 4 quadrants meet in 3x3 region
      char q = 0;
      for (int dj = -1; dj <= 1; ++dj)
      {
        int jdj = j + dj;
        for (int di = -1; di <= 1; ++di)
        {
          int idi = i + di;
          int kdk = idi + jdj * m;
          q |= qs[kdk];
        }
      }
      if (q == 0xF)
      {
        // 4 quadrants meet, check for local minimum at center
        double minmz = 1.0/0.0;
        for (int dj = -1; dj <= 1; ++dj)
        {
          int jdj = j + dj;
          for (int di = -1; di <= 1; ++di)
          {
            int idi = i + di;
            int kdk = idi + jdj * m;
            double mz = cabs(zs[kdk]);
            minmz = mz < minmz ? mz : minmz;
          }
        }
        int k = i + j * m;
        double mz = cabs(zs[k]);
        if (mz <= minmz && minmz < 1.0/0.0)
        {
          // we found a probable zero, output it
          double _Complex ic = ics[k];
          int out;
          #pragma omp atomic capture
          out = o++;
          ops[out] = ip;
          ocs[out] = ic;
        }
      }
    }
  }
  return o;
}

#ifdef IMAGE_E
bool accept_period(m_period_filter_t *data, int p)
{
  (void) data;
  return p >= 129 && (p % 4) != 1;
}

bool reject_period(m_period_filter_t *data, int p)
{
  (void) data;
  return p < 129;
}
#endif

#ifdef IMAGE_F
bool accept_period(m_period_filter_t *data, int p)
{
  (void) data;
  return p >  129 && (p % 4) == 2;
}

bool reject_period(m_period_filter_t *data, int p)
{
  (void) data;
  return p <= 129 || (p % 4) == 1;
}
#endif

int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  const char *filename = "o.png";

  double _Complex center = -1.9409856638234979201929821105592e+00 + I * 6.4820396157412255547279397439918e-04;
#ifdef IMAGE_F
  double radius = 5e-11/2;
#else
  double radius = 5e-9;
#endif

  int maxp = 2000;
  
  int maxiters = 10000;

  int m = 1920;
  int n = 1080;
  int mn = m * n;

  double minfontsize = 1;
  double maxfontsize = 64;

  double er = 25;

  // colouring
  double phi = (sqrt(5) + 1) / 2;
  double gold = 1 / (phi * phi);
  m_pixel_t white = m_pixel_rgba(1, 1, 1, 1);

  // allocate buffers
  double _Complex *cs = malloc(mn * sizeof(*cs));
  double _Complex *zs = malloc(mn * sizeof(*zs));
  char *qs = malloc(mn * sizeof(*qs));
  double _Complex *ocs = malloc(2 * mn * sizeof(*ocs));
  int *ops = malloc(2 * mn * sizeof(*ops));
  m_image *mask = m_image_new(m, n);
  m_image *img = m_image_new(m, n);

  // set image parameters
  m_period_filter_t filter = { &accept_period, &reject_period };
  m_d_transform *t = m_d_transform_rectangular(m, n, center, radius);
  m_d_colour_t *colour = m_d_colour_domain(gold / 5, 0.75, 0.75, &filter);

  // periodicity scan
  initialize_cs(m, n, t, cs);
  memset(zs, 0, mn * sizeof(*zs));
  int o = 0;
  for (int p = 1; p < maxp && o < mn; ++p)
  {
    step_zs(mn, qs, zs, cs);
    o += scan_for_zeroes(m, n, p, ops + o, ocs + o, qs, zs, cs);
  }

  // render fractal
  m_d_render_scanline_filtered(img, t, er, maxiters, colour, &filter);
  m_image_dirty(img);

  // clear mask
  cairo_t *maskcr = cairo_create(m_image_surface(mask));
  cairo_set_source_rgba(maskcr, 1, 1, 1, 1);
  cairo_paint(maskcr);
  cairo_set_source_rgba(maskcr, 0, 0, 0, 1);

  // setup for labels
  cairo_t *cr = cairo_create(m_image_surface(img));
  cairo_set_source_rgba(cr, 1, 1, 1, 1);
  cairo_select_font_face(cr, "LMSans10", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

  // draw labels
  for (int i = 0; i < o; ++i)
  {
    // atom location to image coordinates
    int p = ops[i];
    double _Complex c = ocs[i];
    double _Complex px = c;
    double _Complex dpx = 1;
    m_d_transform_reverse(t, &px, &dpx);
    double x = creal(px);
    double y = cimag(px);
    if (m_image_in_bounds(mask, x, y) && m_image_peek(mask, x, y) == white)
    {
      // no atom already at this location
      // we are the lowest period so far
      // so not a multiple, have true period
      if (m_converged == m_d_nucleus(&c, c, p, 16))
      {
        // check again, in case nucleus is far from located pixel
        px = c;
        if (accept_period(0, p))
        {
          dpx = m_d_filtered_domain_size(c, p, &filter);
//fprintf(stderr, "%d %e\n", p, creal(dpx));
          if (! (creal(dpx) < 1.0/0.0))
            dpx = 0;
          dpx *= 0.5;
        }
        else
        {
          dpx = m_d_domain_size(c, p);
        }
        m_d_transform_reverse(t, &px, &dpx);
        x = creal(px);
        y = cimag(px);
        if (m_image_in_bounds(mask, x, y) && m_image_peek(mask, x, y) == white)
        {
          // calculate font size based on atom domain size
          m_shape shape = m_d_shape_discriminant(m_d_shape_estimate(c, p));
          double fs = (shape == m_cardioid ? 1 : 0.5) * cabs(dpx);
          if (p == 1)
            fs = maxfontsize * 2; // period 1 domain is infinite
          double fsr = 1;
          if (0)//(p % periodmod) != periodneq)
          {
            fs = maxfontsize + 3 * log2(fs);
            fsr = 1.5;
          }
          fs = fmax(fmin(fs, maxfontsize), minfontsize);
          if (fs > minfontsize)
          {
            cairo_set_font_size(cr, fs);
            // draw text centered on point
            char sp[100];
            snprintf(sp, 100, "%d", p);
            cairo_text_extents_t e;
            cairo_text_extents(cr, sp, &e);
            double tx = x - e.x_bearing - e.width / 2.0;
            double ty = y - e.y_bearing - e.height / 2.0;
            cairo_move_to(cr, tx, ty);
            cairo_show_text(cr, sp);
            cairo_fill(cr);
          }
          // update mask
          cairo_arc(maskcr, x, y, fmax(1, fsr * fs), 0, 6.283185307179586);
          cairo_close_path(maskcr);
          cairo_fill(maskcr);
          m_image_flush(mask);
        }
      }
    }
  }

  // save image
  m_image_save_png(img, filename);

  // cleanup
  m_image_delete(img);
  m_image_delete(mask);
  m_d_transform_delete(t);
  m_d_colour_delete(colour);
  free(cs);
  free(zs);
  free(qs);
  free(ocs);
  free(ops);

  // exit
  return 0;
}
