/*

f_c^p(w) = w
f_c^p'(w) = e^(i t)

wn = f_^q(w)

|Dwn| known by output pixel spacing
find |Dt|
binary search between previous Dt/4 and Dt*4?

*/

#define _POSIX_C_SOURCE 2

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include <mandelbrot-numerics.h>
#include <mandelbrot-symbolics.h>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#if 0
typedef int32_t Z;
#define PRZ "d"
typedef uint32_t N;
#define PRN "u"
#else
typedef int64_t Z;
#define PRZ "ld"
typedef uint64_t N;
#define PRN "lu"
#endif

static inline
uint64_t xorshift64(uint64_t x)
{
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  return x;
}

N gcd(N a, N b)
{
  if (a == b) return a;
  if (a < b) return gcd(b, a);
  if (b == 0) return a;
  return gcd(b, a % b);
}

typedef double R;

void hsv2rgb(R h, R s, R v, R *rp, R *gp, R *bp) {
  if (isnan(h) || isinf(h))
  {
    h = 0;
    s = 0;
    v = 0;
  }
  R i, f, p, q, t, r, g, b;
  Z ii;
  if (s == 0.0) {
    r = v;
    g = v;
    b = v;
  } else {
    h = h - floor(h);
    h = h * 6;
    i = floor(h);
    ii = i;
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
      case 0:  r = v; g = t; b = p; break;
      case 1:  r = q; g = v; b = p; break;
      case 2:  r = p; g = v; b = t; break;
      case 3:  r = p; g = q; b = v; break;
      case 4:  r = t; g = p; b = v; break;
      default: r = v; g = p; b = q; break;
    }
  }
  *rp = r;
  *gp = g;
  *bp = b;
}

typedef double _Complex C;

typedef void (*output_t)(void *arg, C w, C c, R t);

#define gain 65536.0
#define line_width 1.0

void trace(C nucleus, N period, N iteration, R Dwn, R scale, output_t output, void *arg)
{
  bool positive = Dwn > 0;
  R t = positive ? 0.5 : -0.5, Dt = positive ? -1e-6 : 1e-6;
  Dt *= scale;
  Dt = fmin(fmax(Dt, -0.5/4), 0.5/4);
  Dwn = fabs(Dwn);
  C w0 = 0, c = nucleus;
  m_d_interior(&w0, &c, w0, c, cexp(2 * M_PI * I * t), period, 64);
  do
  {
    C w00 = w0, c0 = c;
    if (m_failed == m_d_interior(&w0, &c, w0, c, cexp(2 * M_PI * I * t), period, 64))
    {
      break;
    }
    C wn = w0, wn0 = w00;
    for (N i = 0; i < iteration; ++i)
    {
      wn0 = wn0 * wn0 + c0;
      wn = wn * wn + c;
    }
    if (cabs(wn) > 2)
    {
      break;
    }
    if (cabs(wn0 - wn) > 4 * Dwn) break;
    output(arg, wn, c, t);
    R t_min = t;
    R t_max = t + Dt * 4;
    for (N step = 0; step < 16; ++step)
    {
      R t_mid = (t_min + t_max) / 2;
      C w0_next = w0, c_next = c;
      if (m_failed == m_d_interior(&w0_next, &c_next, w0_next, c_next, cexp(2 * M_PI * I * t_mid), period, 64))
      {
        t_max = t_mid;
        continue;
      }
      C wn_next = w0_next;
      for (N i = 0; i < iteration; ++i)
      {
        wn_next = wn_next * wn_next + c;
      }
      R Dwn_next = cabs(wn_next - wn);
      if (Dwn_next > Dwn)
      {
        t_max = t_mid;
      }
      else
      {
        t_min = t_mid;
      }
    }
    R t_next = (t_min + t_max) / 2;
    Dt = t_next - t;
    Dt = fmin(fmax(Dt, -0.5/4), 0.5/4);
    t += Dt;
  }
  while (((positive && t > 0) || (! positive && t < 0)) && fabs(Dt) > 1e-12);
}

struct path
{
  C scale;
  C offset;
  C w, c, t;
  N period;
  bool first;
  bool last;
  N width;
  N height;
  N *histogram;
  N npalette;
  R *palette;
  N aa;
  N prng;
  R r, g, b;
  N depth;
};

void output(void *arg, C w, C c, R t)
{
  struct path *p = arg;
  C nw = p->scale * w + p->offset;
  C pw = p->scale * p->w + p->offset;
  if (! p->first && (p->last || cabs(nw - pw) < 4))
  {
    R desired_points = p->aa * cabs(nw - pw);
    N actual_points = ceil(desired_points);
    R weight = cabs(p->c - c) / cabs(p->w - w); // by length change
    weight *= cabs(p->c - c) * p->scale; // by c length
#if 0
    weight *= weight; // by area
#endif
    weight /= p->period; // by period
    weight /= actual_points; // by points
    weight *= desired_points / actual_points;
    // integrate palette via summed-area table
    Z c0 = floor(fmin(p->t, t) * p->npalette);
    Z c1 = floor(fmax(p->t, t) * p->npalette);
    if (c0 < 0)
    {
      c0 += p->npalette;
      c1 += p->npalette;
    }
    c1 += 1;
    assert(0 <= c0);
    assert(c0 < (Z) p->npalette + 1);
    assert(0 < c1);
    assert(c1 <= (Z) p->npalette + 1);
    R r = (p->palette[3*c1+0] - p->palette[3*c0+0]) / (c1 - c0);
    R g = (p->palette[3*c1+1] - p->palette[3*c0+1]) / (c1 - c0);
    R b = (p->palette[3*c1+2] - p->palette[3*c0+2]) / (c1 - c0);
//    hsv2rgb((p->t + t) * 0.5, 1, 1, &r, &g, &b);
    r += p->r;
    g += p->g;
    b += p->b;
    N depth = p->depth + 1;
    N a  = ceil(gain * weight);
    N ri = ceil(gain * weight * r / depth);
    N gi = ceil(gain * weight * g / depth);
    N bi = ceil(gain * weight * b / depth);
    for (N i = 0; i < actual_points; ++i)
    {
      R l = (i + 0.5) / actual_points;
      C q = pw + l * (nw - pw);
      Z x = creal(q) + line_width * ((p->prng = xorshift64(p->prng)) / 18446744073709551616.0 - 0.5);
      Z y = cimag(q) + line_width * ((p->prng = xorshift64(p->prng)) / 18446744073709551616.0 - 0.5);
      if (0 <= x && x < (Z)p->width && 0 <= y && y < (Z)p->height)
      {
        N k = (y * p->width + x) * 4;
        #pragma omp atomic
        p->histogram[k+0] += a;
        #pragma omp atomic
        p->histogram[k+1] += ri;
        #pragma omp atomic
        p->histogram[k+2] += gi;
        #pragma omp atomic
        p->histogram[k+3] += bi;
      }
    }
  }
  p->first = false;
  p->w = w;
  p->c = c;
  p->t = t;
}

N trace_tree(struct path *p, C rootw, C rootc, C nucleus, N period, N max_period, N depth, char *provenance, char *provenance_end)
{
  N total = 0;
  R size = cabs(m_d_size(nucleus, period));
  for (N branch = 0; branch < period; ++branch)
  {
    if (false)
    {
      fprintf(stderr, "%14"PRN " : %s : %"PRN "\n", total, provenance, branch);
    }
    R r0 = p->r, g0 = p->g, b0 = p->b, r1, g1, b1;
    if (period > 1)
    {
      hsv2rgb(branch / (R) period, 1, 1, &r1, &g1, &b1);
      p->r = r0 + r1;
      p->g = g0 + g1;
      p->b = b0 + b1;
      p->depth += 1;
    }
    p->period = period;
    p->first = true;
    p->last = false;
    trace(nucleus, period, branch, +1 / cabs(p->scale), 1 / size, output, p);
    p->last = true;
    output(p, rootw, rootc, 0);
    p->first = true;
    p->last = false;
    trace(nucleus, period, branch, -1 / cabs(p->scale), 1 / size, output, p);
    p->last = true;
    output(p, rootw, rootc, 0);
    if (period > 1)
    {
      p->depth -= 1;
      p->r = r0;
      p->g = g0;
      p->b = b0;
    }
    rootw = rootw * rootw + rootc;
    ++total;
  }
  for (N den = 2; den * period <= max_period; ++den)
  {
    for (N num = 1; num < den; ++num)
    {
      if (gcd(num, den) == 1)
      {
        sprintf(provenance_end, " %" PRN "/%"PRN " %"PRN, num, den, den * period);
        C c = nucleus, w = 0;
        m_d_interior(&w, &c, w, c, cexp(2 * M_PI * I * num / den), period, 64);
        C wp = w, dwp = 1, ddwp = 0;
        for (N i = 0; i < period; ++i)
        {
          ddwp = 2 * (wp * ddwp + dwp * dwp);
          dwp = 2 * wp * dwp;
          wp = wp * wp + c;
        }
        bool cardioid = period == 1;
        C child = c + ddwp / cabs(ddwp) * size / (den * den) * (cardioid ? sin(M_PI * num / den) : 1);
        m_d_nucleus(&child, child, den * period, 64);
        R r0 = p->r, g0 = p->g, b0 = p->b, r1, g1, b1;
        hsv2rgb(num / (R) den, 1, 1, &r1, &g1, &b1);
        p->r = r0 + r1;
        p->g = g0 + g1;
        p->b = b0 + b1;
        p->depth += 1;
        total += trace_tree(p, w, c, child, den * period, max_period, depth + 1, provenance, provenance_end + strlen(provenance_end));
        p->depth -= 1;
        p->r = r0;
        p->g = g0;
        p->b = b0;
      }
    }
  }
  return total;
}

void tonemap(const N *histogram, uint8_t *ppm, N w, N h)
{
  N total = 0;
  for (N y = 0; y < h; ++y)
  {
    for (N x = 0; x < w; ++x)
    {
      N k = (y * w + x) * 4;
      total += histogram[k+0];
    }
  }
  R s = w * h / (R) total;
  #pragma omp parallel for
  for (N y = 0; y < h; ++y)
  {
    for (N x = 0; x < w; ++x)
    {
      N k = (y * w + x) * 4;
      N j = (y * w + x) * 3;
#if 1
      R a = histogram[k + 0];
      R r = 0, g = 0, b = 0;
      if (a > 0)
      {
        R v = log(1 + log(1 + s * a)) / a;
        r = 195 * sqrt(v * histogram[k+1]);
        g = 195 * sqrt(v * histogram[k+2]);
        b = 195 * sqrt(v * histogram[k+3]);
      }
#else
      R r = 256 * sqrt(sqrt(histogram[k+1] / gain));
      R g = 256 * sqrt(sqrt(histogram[k+2] / gain));
      R b = 256 * sqrt(sqrt(histogram[k+3] / gain));
#endif
      ppm[j+0] = fmin(fmax(r, 0), 255);
      ppm[j+1] = fmin(fmax(g, 0), 255);
      ppm[j+2] = fmin(fmax(b, 0), 255);
    }
  }
}

int main(int argc, char **argv)
{
  N max_period = argc > 1 ? atoi(argv[1]) : 3;
  N sharpness = 8;
#if 0
  N w = 4096;
  N h = 4096;
  R r = 2;
#else
  N w = 1920;
  N h = 1080;
  R r = 1.25;
#endif
  N *histogram = malloc(sizeof(*histogram) * w * h * 4);
  N npalette = 1024;
  R *palette = malloc(sizeof(*palette) * (npalette + 2) * 3);
  palette[0] = 0;
  palette[1] = 0;
  palette[2] = 0;
  for (N i = 1; i <= npalette + 1; ++i)
  {
    hsv2rgb((i - 0.5) / npalette, 1, 1, &palette[3*i+0], &palette[3*i+1], &palette[3*i+2]);
    palette[3*i+0] += palette[3*(i-1)+0];
    palette[3*i+1] += palette[3*(i-1)+1];
    palette[3*i+2] += palette[3*(i-1)+2];
  }
  N total = 0;
  N jobs = 0;
  for (N period = 1; period <= max_period; ++period)
  {
    jobs += ((N) 1 << (period - 1));
  }
  N progress = 0;
  #pragma omp parallel for schedule(dynamic, 1)
  for (N job = 0; job < jobs; ++job)
  {
    bool have_job = false;
    N period, num, den;
    N myjob = job;
    for (N p = 1; p <= max_period; ++p)
    {
      if (myjob < (N) 1 << (p - 1))
      {
        period = p;
        num = (myjob << 1) + 1;
        den = (1 << p) - 1;
        have_job = true;
        break;
      }
      else
      {
        myjob -= 1 << (p - 1);
      }
    }
    if (have_job)
    {
      if (num >= den) num -= den;
      N prng = job + 1;
      // FIXME check if angle corresponds to minibrot or child
      mpq_t e;
      mpq_init(e);
      mpq_set_ui(e, num, den);
      mpq_canonicalize(e);
      if (m_period(e) == (Z) period)
      {
        C nucleus = m_d_exray_in_do(e, sharpness, 4 * sharpness * period, 64);
        mpq_clear(e);
        m_d_nucleus(&nucleus, nucleus, period, 64);
        if (m_cardioid == m_d_shape(nucleus, period))
        {
          C rootw = 0, rootc = nucleus;
          m_d_interior(&rootw, &rootc, rootw, rootc, 1, period, 64);
          struct path p =
            { h / (2 * r), (w + h * I) / 2.0
            , 0, 0, 0, 0, 0, 0
            , w, h, histogram
            , npalette, palette
            , 64, prng
            , 0, 0, 0, 0
            };
          char provenance[1000];
          sprintf(provenance, "%" PRN, period);
          N my_total = trace_tree(&p, rootw, rootc, nucleus, period, max_period, 0, provenance, provenance + strlen(provenance));
          prng = p.prng;
          #pragma omp atomic
          total += my_total;
        }
      }
    }
    N myprogress;
    #pragma omp atomic capture
    myprogress = ++progress;
    if ((myprogress - 1) * 100 / jobs != myprogress * 100 / jobs)
    {
      #pragma omp critical
      fprintf(stderr, "%4"PRN"%%\r", myprogress * 100 / jobs);
    }
  }
  fprintf(stderr, "\n%"PRN"\n", total);
  uint8_t *ppm = malloc(sizeof(*ppm) * w * h * 3);
  tonemap(histogram, ppm, w, h);
  printf("P6\n%"PRN " %"PRN "\n255\n", w, h);
  fwrite(ppm, sizeof(*ppm) * w * h * 3, 1, stdout);
  return 0;
}
