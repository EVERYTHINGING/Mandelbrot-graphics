#include <mandelbrot-graphics.h>
#include "m_d_util.h"

extern void m_d_render_scanline_filtered(m_image *image, m_d_transform *transform, double er, int maxiters, m_d_colour_t *colour, m_period_filter_t *filter)
{
  m_image_flush(image);
  int width = m_image_get_width(image);
  int height = m_image_get_height(image);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < height; ++j)
  {
    m_d_compute *px = m_d_compute_alloc(maxiters);
    for (int i = 0; i < width; ++i)
    {
      double _Complex c = i + I * j;
      double _Complex dc = 1;
      m_d_transform_forward(transform, &c, &dc);
      m_d_compute_init(px, m_d_compute_get_tag(px), er, c, filter);
      m_d_compute_step(px, maxiters);
      m_image_plot(image, i, j, m_d_colour(colour, px, er, dc));
    }
    m_d_compute_free(px);
  }
  m_image_dirty(image);
}

extern void m_d_render_scanline(m_image *image, m_d_transform *transform, double er, int maxiters, m_d_colour_t *colour)
{
  m_d_render_scanline_filtered(image, transform, er, maxiters, colour, 0);
}
