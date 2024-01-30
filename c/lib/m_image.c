#include <mandelbrot-graphics.h>
#include <stdio.h>
#include "m_d_util.h"

struct m_image {
  int width;
  int height;
  int stride;
  m_pixel_t *data;
  cairo_surface_t *surface;
};

extern m_image *m_image_new(int width, int height) {
  m_image *img = (m_image *) calloc(1, sizeof(*img));
  img->width = width;
  img->height = height;
  img->stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, width);
  img->data = (m_pixel_t *) calloc(1, img->stride * height);
  img->surface = cairo_image_surface_create_for_data((unsigned char *) img->data, CAIRO_FORMAT_ARGB32, width, height, img->stride);
  return img;
}

extern void m_image_delete(m_image *img) {
  cairo_surface_destroy(img->surface);
  free(img->data);
  free(img);
}

extern m_image *m_image_load_ppm(const char *filename) {
  FILE *in = fopen(filename, "rb");
  int w = 0, h = 0;
  fscanf(in, "P6\n%d %d 255", &w, &h);
  if (! (w > 0 && h > 0 && '\n' == fgetc(in))) {
    fclose(in);
    return 0;
  }
  unsigned char *ppm = malloc(3 * w * h);
  if (! (ppm && 1 == fread(ppm, 3 * w * h, 1, in))) {
    fclose(in);
    return 0;
  }
  fclose(in);
  m_image *img = m_image_new(w, h);
  if (! img) {
    free(ppm);
    return 0;
  }
  #pragma omp parallel for
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      m_image_plot(img, i, j, m_pixel_rgba
        ( ppm[3 * (w * j + i) + 0] / 255.0
        , ppm[3 * (w * j + i) + 1] / 255.0
        , ppm[3 * (w * j + i) + 2] / 255.0
        , 1.0
        ));
    }
  }
  m_image_dirty(img);
  free(ppm);
  return img;
}

extern cairo_surface_t *m_image_surface(m_image *img) {
  return img->surface;
}

extern void m_image_flush(m_image *img) {
  cairo_surface_flush(img->surface);
}

extern void m_image_dirty(m_image *img) {
  cairo_surface_mark_dirty(img->surface);
}

extern bool m_image_in_bounds(m_image *img, int x, int y) {
  return 0 <= x && x < img->width && 0 <= y && y < img->height;
}

extern void m_image_plot(m_image *img, int x, int y, m_pixel_t c) {
  img->data[y * (img->stride >> 2) + x] = c;
}

extern m_pixel_t m_image_peek(m_image *img, int x, int y) {
  return img->data[y * (img->stride >> 2) + x];
}

extern bool m_image_save_png(m_image *img, const char *filename) {
  return CAIRO_STATUS_SUCCESS == cairo_surface_write_to_png(img->surface, filename);
}

extern int m_image_get_width(const m_image *img) {
  return img->width;
}

extern int m_image_get_height(const m_image *img) {
  return img->height;
}
