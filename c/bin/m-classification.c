#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mandelbrot-graphics.h>

#define UNKNOWN 0
#define INTERIOR 1
#define EXTERIOR 2

int classify(double d, double s)
{
#define SQRT2 1.4142135623730951
	int k = 0;
	if (d < 0)
	{
		k |= (1<<INTERIOR);
		if (-s < 4 * d)
			k |= (1<<EXTERIOR);
		else if (-4 * SQRT2 * s < d)
			k |= (1<<UNKNOWN);
	}
	else
	{
		k |= (1<<EXTERIOR);
		if (4 * d < s)
			k |= (1<<INTERIOR);
		else if (d < 4 * SQRT2 * s)
			k |= (1<<UNKNOWN);
	}
	return k;
#undef SQRT2
}

#ifdef USE_ASCII
#define HEIGHT 50
#define WIDTH 100
#define PAR 0.5
#else
#define HEIGHT 1080
#define WIDTH 1920
#define PAR 1
#endif

char palette[2][2][2] = {{{' ','X'},{'I','@'}},{{'?','x'},{'i','#'}}};
char image[HEIGHT+1][WIDTH+1];
unsigned char colours[2][2][2][3] =
{{{{128,128,128}
  ,{128,255,255}}
 ,{{255,255,128}
  ,{0,0,0}}}
,{{{255,0,0}
  ,{0,128,128}}
 ,{{128,128,0}
  ,{0,0,0}}}};
unsigned char rgb[HEIGHT][WIDTH][3];

#define ITERATIONS (1<<24)
#define PARTIALS ITERATIONS
#define ER (1e10)

int main(int argc, char **argv)
{
	(void) argc;
	(void) argv;
	double _Complex c0 = -0.75;
	double r = 1.5;
	memset(image, ' ', sizeof(image));
	m_d_compute *m = m_d_compute_alloc(PARTIALS);
	m_compute_t bias = m_exterior;
	for (int SUPER = 1; SUPER <= 256; SUPER <<= 1)
	{
		double s = r * 2 / (SUPER * HEIGHT);
		for (int y = 0; y < HEIGHT; ++y)
		{
			for (int x = 0; x < WIDTH; ++x)
			{
				bias = image[y][x] == 'i' || image[y][x] == 'I' ? m_interior : image[y][x] == 'x' || image[y][x] == 'X' ? m_exterior : bias;
				if (image[y][x] == ' ' || image[y][x] == 'i' || image[y][x] == 'x')
				{
					int p[3] = { 0, 0, 0 };
					for (int j = 0; j < SUPER; ++j)
					{
						for (int i = 0; i < SUPER; ++i)
						{
							double u = ((SUPER * x + i + 0.5) / (SUPER * WIDTH) - 0.5) * 2 * WIDTH / HEIGHT * PAR;
							double v = (0.5 - (SUPER * y + j + 0.5) / (SUPER * HEIGHT)) * 2;
							double _Complex c = c0 + r * (u + I * v);;
							m_d_compute_init(m, bias, ER, c, 0);
							m_d_compute_step(m, ITERATIONS);
							bias = m_d_compute_get_tag(m);
							double d = m_d_compute_get_de(m) * (bias == m_interior ? -1 : +1);
							int k = classify(d, s);
							p[UNKNOWN ] += !!(k&(1<<UNKNOWN ));
							p[INTERIOR] += !!(k&(1<<INTERIOR));
							p[EXTERIOR] += !!(k&(1<<EXTERIOR));
							if (p[INTERIOR] && p[EXTERIOR])
								break;
						}
						if (p[INTERIOR] && p[EXTERIOR])
							break;
					}
					image[y][x] = palette[!!p[UNKNOWN]][!!p[INTERIOR]][!!p[EXTERIOR]];
					for (int i = 0; i < 3; ++i)
						rgb[y][x][i] = colours[!!p[UNKNOWN]][!!p[INTERIOR]][!!p[EXTERIOR]][i];
				}
			}
			image[y][WIDTH] = '\n';
		}
		image[HEIGHT][0] = '\n';
		image[HEIGHT][1] = 0;
#ifdef USE_ASCII
		puts(image);
#else
		printf("P6\n%d %d\n255\n", WIDTH, HEIGHT);
		fwrite(rgb, sizeof(rgb), 1, stdout);
#endif
		fflush(stdout);
	}
	m_d_compute_free(m);
	return 0;
}
