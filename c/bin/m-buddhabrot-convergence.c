#include <math.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <mandelbrot-graphics.h>

bool running = false;

#define class_unknown 0
#define class_exterior 1
#define class_interior 2

struct context {
  uint64_t maxiters;
  uint64_t height;
  uint64_t width;

  double dx;
  double dy;

  double exterior_de_threshold_hi;
  double exterior_de_threshold_lo;
  double exterior_de_decrement;

  double interior_de_threshold_hi;
  double interior_de_threshold_lo;
  double interior_de_decrement;

  uint64_t last_class;
  double last_exterior_de;
  double last_interior_de;
  uint64_t last_interior_period;

  uint64_t i;
  uint64_t j;

  uint64_t exterior_without_de;
  uint64_t exterior_with_de;
  uint64_t boundary_with_interior_and_de;
  uint64_t interior_with_de;
  uint64_t interior_without_de;

  uint64_t interior;
  uint64_t unknown;
  uint64_t level_set_count[];
};

void get_c(struct context *ctx, double *cx, double *cy) {
  *cx = (ctx->i + ctx->dx) / ctx->width * 4.04 - 2.02;
  *cy = (ctx->j + ctx->dy) / ctx->height * 2.02;
}

void advance_ij(struct context *ctx) {
  ++(ctx->i);
  if (ctx->i == ctx->width) {
    ctx->i = 0;
    ++(ctx->j);
    ctx->last_class = class_unknown;
  }
}

void compute_exterior_pixel_without_de(struct context *ctx) {
  const double er2 = 4;
  uint64_t maxiters = ctx->maxiters;
  double cx, cy;
  get_c(ctx, &cx, &cy);
  double zx = cx, zy = cy, zx2, zy2, z2, zxy;
  uint64_t n = 0;
  do {
    zx2 = zx * zx;
    zy2 = zy * zy;
    z2 = zx2 + zy2;
    if (z2 > er2) { break; }
    zxy = zx * zy;
    zx = zx2 - zy2 + cx;
    zy = 2 * zxy + cy;
    ++n;
  } while (n < maxiters);
  if (n < maxiters) {
    ++(ctx->level_set_count[n]);
    ctx->last_class = class_exterior;
    ctx->last_exterior_de -= ctx->exterior_de_decrement;
  } else {
    ++(ctx->unknown);
    ctx->last_class = class_unknown;
  }
  ++(ctx->exterior_without_de);
}

void compute_exterior_pixel_with_de(struct context *ctx) {
  const double er2_lo = 4;
  const double er2_hi = 65536;
  uint64_t maxiters = ctx->maxiters;
  double cx, cy;
  get_c(ctx, &cx, &cy);
  double zx = cx, zy = cy, zx2, zy2, z2, zxy, dzx = 1, dzy = 0, dzt;
  uint64_t n = 0;
  uint64_t m = maxiters;
  do {
    zx2 = zx * zx;
    zy2 = zy * zy;
    z2 = zx2 + zy2;
    if (z2 > er2_hi) {
      break;
    }
    if (z2 < er2_lo) {
      m = n;
    }
    dzt = 2 * (zx * dzx - zy * dzy) + 1;
    dzy = 2 * (zx * dzy + zy * dzx);
    dzx = dzt;
    zxy = zx * zy;
    zx = zx2 - zy2 + cx;
    zy = 2 * zxy + cy;
    ++n;
  } while (n < maxiters);
  ++m;
  if (n < maxiters && m < maxiters) {
    ++(ctx->level_set_count[m]);
    ctx->last_class = class_exterior;
    ctx->last_exterior_de = sqrt(z2) * log(z2) / hypot(dzx, dzy);
  } else {
    ++(ctx->unknown);
    ctx->last_class = class_unknown;
  }
  ++(ctx->exterior_with_de);
}

void compute_boundary_pixel_with_interior_checking_and_de(struct context *ctx) {
  const double er2_lo = 4;
  const double er2_hi = 65536;
  uint64_t maxiters = ctx->maxiters;
  double cx, cy;
  get_c(ctx, &cx, &cy);
  double zx = cx, zy = cy, zx2, zy2, z2, zxy, dzx = 1, dzy = 0, dzt, mz2 = 1.0 / 0.0;
  double interior_de = -1;
  uint64_t n = 0;
  uint64_t m = maxiters;
  uint64_t result = class_unknown;
  uint64_t period = 0;
  do {
    zx2 = zx * zx;
    zy2 = zy * zy;
    z2 = zx2 + zy2;
    if (z2 > er2_hi) {
      result = class_exterior;
      break;
    }
    if (z2 < er2_lo) {
      m = n;
    }
    if (z2 < mz2) {
      mz2 = z2;
      interior_de = -1;
      complex double dz_unused = 0;
      if (m_d_interior_de(&interior_de, &dz_unused, zx + I * zy, cx + I * cy, n + 1, 64)) {
        result = class_interior;
        period = n + 1;
        break;
      }
    }
    dzt = 2 * (zx * dzx - zy * dzy) + 1;
    dzy = 2 * (zx * dzy + zy * dzx);
    dzx = dzt;
    zxy = zx * zy;
    zx = zx2 - zy2 + cx;
    zy = 2 * zxy + cy;
    ++n;
  } while (n < maxiters);
  ++m;
  if (! (m < maxiters)) {
    result = class_unknown;
  }
  switch (result) {
  case class_interior:
    ++(ctx->interior);
    ctx->last_interior_de = interior_de;
    ctx->last_interior_period = period;
    break;
  case class_exterior:
    ++(ctx->level_set_count[m]);
    ctx->last_class = class_exterior;
    ctx->last_exterior_de = sqrt(z2) * log(z2) / hypot(dzx, dzy);
    break;
  default:
    ++(ctx->unknown);
    ctx->last_class = class_unknown;
    break;
  }
  ++(ctx->boundary_with_interior_and_de);
}

void compute_interior_pixel_without_de(struct context *ctx) {
  ++(ctx->interior);
  ctx->last_interior_de -= ctx->interior_de_decrement;
  ++(ctx->interior_without_de);
}

void compute_interior_pixel_with_de(struct context *ctx) {
  double cx, cy;
  get_c(ctx, &cx, &cy);
  complex double c = cx + I * cy;
  complex double z = 0;
  complex double dz_unused;
  uint64_t period = ctx->last_interior_period;
  for (uint64_t i = 0; i < period; ++i) {
    z = z * z + c;
  }
  double interior_de = -1;
  if (m_d_interior_de(&interior_de, &dz_unused, z, c, period, 64)) {
    ++(ctx->interior);
    ctx->last_class = class_interior;
    ctx->last_interior_de = interior_de;
  } else {
    ++(ctx->unknown);
    ctx->last_class = class_unknown;
  }
  ++(ctx->interior_with_de);
}

int main(int argc, char **argv) {
  if (argc > 1) {
    if (strcmp("init", argv[1]) == 0) {
      if (! (argc > 4)) { return 1; }
      uint64_t maxiters = ((uint64_t) 1) << atoi(argv[2]);
      uint64_t height = ((uint64_t) 1) << atoi(argv[3]);
      uint64_t width = 2 * height;
      double pixel_size = 2.02 / height;

      struct context *ctx = calloc(1, sizeof(*ctx));
      ctx->maxiters = maxiters;
      ctx->height = height;
      ctx->width = width;
      ctx->dx = 0.5;
      ctx->dy = 0.5;
      ctx->exterior_de_threshold_hi = 8 * pixel_size;
      ctx->exterior_de_threshold_lo = 4 * pixel_size;
      ctx->exterior_de_decrement = pixel_size;
      ctx->interior_de_threshold_hi = 8 * pixel_size;
      ctx->interior_de_threshold_lo = 4 * pixel_size;
      ctx->interior_de_decrement = pixel_size;

      FILE *out = fopen(argv[4], "wb");
      if (! out) { return 1; }
      if (1 != fwrite(ctx, sizeof(*ctx), 1, out)) { return 1; }
      uint64_t bytes = sizeof(struct context) + maxiters * sizeof(uint64_t);
      uint64_t zero = 0;
      if (-1 == fseek(out, bytes - sizeof(zero), SEEK_SET)) { return 1; }
      if (1 != fwrite(&zero, sizeof(zero), 1, out)) { return 1; }
      fclose(out);
      return 0;
    }

    if (! (argc > 2)) { return 1; }
    int fd = open(argv[2], O_RDWR);
    uint64_t maxiters = 0;
    read(fd, &maxiters, sizeof(maxiters));
    lseek(fd, 0, SEEK_SET);
    uint64_t bytes = sizeof(struct context) + maxiters * sizeof(uint64_t);
    struct context *ctx = mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    if (strcmp("result", argv[1]) == 0) {
      if (ctx->j < ctx->height) {
        fprintf(stderr, "INCOMPLETE\n");
      } else {
        printf("# levelset area\n");
        for (uint64_t m = 0; m < ctx->maxiters; ++m) {
          double area = ctx->level_set_count[m] * 4.04 / ctx->width * 2.02 / ctx->height;
          if (area > 0) {
            printf("%lu %.16e\n", m, area);
          }
        }
      }
    } else if (strcmp("compute", argv[1]) == 0) {
      running = true;
      while (running && ctx->j < ctx->height) {
        if (ctx->last_class == class_exterior) {
          if (ctx->last_exterior_de > ctx->exterior_de_threshold_hi) {
            compute_exterior_pixel_without_de(ctx);
          } else if (ctx->last_exterior_de > ctx->exterior_de_threshold_lo) {
            compute_exterior_pixel_with_de(ctx);
          } else {
            compute_boundary_pixel_with_interior_checking_and_de(ctx);
          }
        } else if (ctx->last_class == class_interior) {
          if (ctx->last_interior_de > ctx->interior_de_threshold_hi) {
            compute_interior_pixel_without_de(ctx);
          } else if (ctx->last_interior_de > ctx->interior_de_threshold_lo) {
            compute_interior_pixel_with_de(ctx);
          } else {
            compute_boundary_pixel_with_interior_checking_and_de(ctx);
          }
        } else {
          compute_boundary_pixel_with_interior_checking_and_de(ctx);
        }
        advance_ij(ctx);
      }
      if (! (ctx->j < ctx->height)) {
        printf("DONE\n");
      }
    }

    munmap(ctx, bytes);
    return 0;
  }
}
