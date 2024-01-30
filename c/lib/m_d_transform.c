#include <mandelbrot-graphics.h>
#include "m_d_util.h"

extern void m_d_transform_forward(m_d_transform *t, double _Complex *c, double _Complex *dc) {
  t->forward(t, c, dc);
}

extern void m_d_transform_reverse(m_d_transform *t, double _Complex *c, double _Complex *dc) {
  t->reverse(t, c, dc);
}

extern void m_d_transform_delete(m_d_transform *t) {
  free(t);
}

struct m_d_transform_compose_t {
  m_d_transform transform;
  m_d_transform *s;
  m_d_transform *t;
};
typedef struct m_d_transform_compose_t m_d_transform_compose_t;

static void m_d_transform_compose_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_compose_t *t = userdata;
  m_d_transform_forward(t->s, c, dc);
  m_d_transform_forward(t->t, c, dc);
}

static void m_d_transform_compose_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_compose_t *t = userdata;
  m_d_transform_reverse(t->t, c, dc);
  m_d_transform_reverse(t->s, c, dc);
}

extern m_d_transform *m_d_transform_compose(m_d_transform *s, m_d_transform *t) {
  m_d_transform_compose_t *r = malloc(sizeof(*r));
  r->transform.forward = m_d_transform_compose_forward;
  r->transform.reverse = m_d_transform_compose_reverse;
  r->s = s;
  r->t = t;
  return &r->transform;
}

extern m_d_transform *m_d_transform_invert(m_d_transform *t) {
  m_d_transform_f *f = t->forward;
  t->forward = t->reverse;
  t->reverse = f;
  return t;
}

struct m_d_transform_rectangular_t {
  m_d_transform transform;
  int width;
  int height;
  double _Complex c;
  double r;
};
typedef struct m_d_transform_rectangular_t m_d_transform_rectangular_t;

static void m_d_transform_rectangular_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_rectangular_t *rect = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = rect->c + rect->r * conj(c0 - ((rect->width/2 - 0.5) + I * (rect->height/2 - 0.5))) / (rect->height/2);
  *dc = conj(dc0 * rect->r / (rect->height/2));
}

static void m_d_transform_rectangular_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_rectangular_t *rect = (m_d_transform_rectangular_t *) userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = conj((c0 - rect->c) / rect->r * (rect->height/2)) + ((rect->width/2 - 0.5) + I * (rect->height/2 - 0.5));
  *dc = conj(dc0 / rect->r * (rect->height/2));
}

extern m_d_transform *m_d_transform_rectangular(int w, int h, double _Complex c, double r) {
  m_d_transform_rectangular_t *rect = malloc(sizeof(*rect));
  rect->transform.forward = m_d_transform_rectangular_forward;
  rect->transform.reverse = m_d_transform_rectangular_reverse;
  rect->width = w;
  rect->height = h;
  rect->c = c;
  rect->r = r;
  return &rect->transform;
}


struct m_d_transform_linear_t {
  m_d_transform transform;
  double _Complex add;
  double _Complex mul;
};
typedef struct m_d_transform_linear_t m_d_transform_linear_t;

static void m_d_transform_linear_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_linear_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = t->mul * c0 + t->add;
  *dc = t->mul * dc0;
}

static void m_d_transform_linear_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_linear_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = (c0 - t->add) / t->mul;
  *dc = dc0 / t->mul;
}

extern m_d_transform *m_d_transform_linear(double _Complex add, double _Complex mul) {
  m_d_transform_linear_t *t = malloc(sizeof(*t));
  t->transform.forward = m_d_transform_linear_forward;
  t->transform.reverse = m_d_transform_linear_reverse;
  t->add = add;
  t->mul = mul;
  return &t->transform;
}

struct m_d_transform_moebius_t {
  m_d_transform transform;
  m_d_mat2 m;
};
typedef struct m_d_transform_moebius_t m_d_transform_moebius_t;

static void m_d_transform_moebius_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_moebius_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  double _Complex d = t->m.c * c0 + t->m.d;
  *c = (t->m.a * c0 + t->m.b) / d;
  *dc = dc0 * m_d_mat2_det(&t->m) / (d * d);
}

static void m_d_transform_moebius_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_moebius_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  double _Complex d = -t->m.c * c0 + t->m.a;
  *c = (t->m.d * c0 - t->m.b) / d;
  *dc = dc0 * m_d_mat2_det(&t->m) / (d * d);
}

extern m_d_transform *m_d_transform_moebius(const m_d_mat2 *m) {
  m_d_transform_moebius_t *t = malloc(sizeof(*t));
  t->transform.forward = m_d_transform_moebius_forward;
  t->transform.reverse = m_d_transform_moebius_reverse;
  m_d_mat2_set(&t->m, m);
  return &t->transform;
}

extern m_d_transform *m_d_transform_moebius3(double _Complex zero, double _Complex one, double _Complex infinity) {
  m_d_mat2 m;
  m_d_mat2_moebius3(&m, zero, one, infinity);
  return m_d_transform_moebius(&m);
}

extern m_d_transform *m_d_transform_moebius6(double _Complex a0, double _Complex a1, double _Complex a2, double _Complex b0, double _Complex b1, double _Complex b2) {
  m_d_mat2 ma, mb, m;
  m_d_mat2_moebius3(&ma, a0, a1, a2);
  m_d_mat2_moebius3(&mb, b0, b1, b2);
  m_d_mat2_inv(&mb, &mb);
  m_d_mat2_mul(&m, &ma, &mb);
  return m_d_transform_moebius(&m);
}

struct m_d_transform_cardioid_t {
  m_d_transform transform;
};
typedef struct m_d_transform_cardioid_t m_d_transform_cardioid_t;

static void m_d_transform_cardioid_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  (void) userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  double _Complex r = csqrt(1 - 4 * c0);
  *c = r - 1;
  *dc = dc0 * -2 / r;
}

static void m_d_transform_cardioid_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  (void) userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  double _Complex r = c0 + 1;
  *c = (1 - r * r) / 4;
  *dc = dc0 * -r / 2;
}

extern m_d_transform *m_d_transform_cardioid() {
  m_d_transform_cardioid_t *t = malloc(sizeof(*t));
  t->transform.forward = m_d_transform_cardioid_forward;
  t->transform.reverse = m_d_transform_cardioid_reverse;
  return &t->transform;
}

struct m_d_transform_exponential_t {
  m_d_transform transform;
  double _Complex center;
};
typedef struct m_d_transform_exponential_t m_d_transform_exponential_t;

static void m_d_transform_exponential_reverse(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_exponential_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = clog(c0 - t->center);
  *dc = dc0 / (c0 - t->center);
}

static void m_d_transform_exponential_forward(void *userdata, double _Complex *c, double _Complex *dc) {
  m_d_transform_exponential_t *t = userdata;
  double _Complex c0 = *c;
  double _Complex dc0 = *dc;
  *c = cexp(c0) + t->center;
  *dc = dc0 * cexp(c0);
}

extern m_d_transform *m_d_transform_exponential(double _Complex center) {
  m_d_transform_exponential_t *t = malloc(sizeof(*t));
  t->transform.forward = m_d_transform_exponential_forward;
  t->transform.reverse = m_d_transform_exponential_reverse;
  t->center = center;
  return &t->transform;
}
