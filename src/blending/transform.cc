#include <blending/transform.h>

#include <float.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdarg.h>
#include <cstring>
#include <cstdlib>
#include <ctype.h>
#include <limits.h>

namespace blending {
#define R_EPS  1.0e-6
#define MAXITER 100

#ifndef abs
#define abs(a) ( (a) >= 0 ? (a) : -(a) )
#endif

#ifdef _MSC_VER
#define isnan(a) _isnan(a)
#define isinf(a) (_fpclass(a) == _FPCLASS_NINF || _fpclass(a) == _FPCLASS_PINF)
#endif

void    matrix_matrix_mult(double m1[3][3], double m2[3][3], double result[3][3]);
int     polzeros_();

void CubeZero(double *a, int *n, double *root);
void SquareZero(double *a, int *n, double *root);
double CubeRoot(double x);


//------------------------- Some auxilliary math functions --------------------------------------------

void matrix_mult(double m[3][3], double vector[3])
{
	register int i;
	register double v0 = vector[0];
	register double v1 = vector[1];
	register double v2 = vector[2];


	for (i = 0; i<3; i++)
	{
		vector[i] = m[i][0] * v0 + m[i][1] * v1 + m[i][2] * v2;
	}
}

void matrix_inv_mult(double m[3][3], double vector[3])
{
	register int i;
	register double v0 = vector[0];
	register double v1 = vector[1];
	register double v2 = vector[2];

	for (i = 0; i<3; i++)
	{
		vector[i] = m[0][i] * v0 + m[1][i] * v1 + m[2][i] * v2;
	}
}

void matrix_matrix_mult(double m1[3][3], double m2[3][3], double result[3][3])
{
	register int i, k;

	for (i = 0; i<3; i++)
	{
		for (k = 0; k<3; k++)
		{
			result[i][k] = m1[i][0] * m2[0][k] + m1[i][1] * m2[1][k] + m1[i][2] * m2[2][k];
		}
	}
}

// Set matrix elements based on Euler angles a, b, c

void SetMatrix(double a, double b, double c, double m[3][3], int cl)
{
	double mx[3][3], my[3][3], mz[3][3], dummy[3][3];


	// Calculate Matrices;

	mx[0][0] = 1.0;                                mx[0][1] = 0.0;                                mx[0][2] = 0.0;
	mx[1][0] = 0.0;                                mx[1][1] = cos(a);                     mx[1][2] = sin(a);
	mx[2][0] = 0.0;                                mx[2][1] = -mx[1][2];                   mx[2][2] = mx[1][1];

	my[0][0] = cos(b);                              my[0][1] = 0.0;                                my[0][2] = -sin(b);
	my[1][0] = 0.0;                                my[1][1] = 1.0;                                my[1][2] = 0.0;
	my[2][0] = -my[0][2];                   my[2][1] = 0.0;                                my[2][2] = my[0][0];

	mz[0][0] = cos(c);                     mz[0][1] = sin(c);                     mz[0][2] = 0.0;
	mz[1][0] = -mz[0][1];                   mz[1][1] = mz[0][0];                   mz[1][2] = 0.0;
	mz[2][0] = 0.0;                                mz[2][1] = 0.0;                                mz[2][2] = 1.0;

	if (cl)
		matrix_matrix_mult(mz, mx, dummy);
	else
		matrix_matrix_mult(mx, mz, dummy);
	matrix_matrix_mult(dummy, my, m);
}
//------------------------------- Transformation functions --------------------------------------------


#define         distanceparam   (*((double*)params))
#define         shift           (*((double*)params))
#define         var0            ((double*)params)[0]
#define         var1            ((double*)params)[1]
#define         var2            ((double*)params)[2]
#define         var3            ((double*)params)[3]
#define         mp              ((struct TransParams*)params)

// Rotate equirectangular image

int rotate_erect(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double 180degree_turn(screenpoints), double turn(screenpoints);

	*x_src = x_dest + var1;

	while (*x_src < -var0)
		*x_src += 2 * var0;

	while (*x_src >  var0)
		*x_src -= 2 * var0;

	*y_src = y_dest;
	return 1;
}



// Calculate inverse 4th order polynomial correction using Newton
// Don't use on large image (slow)!


int inv_radial(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double coefficients[4]

	register double rs, rd, f, scale;
	int iter = 0;

	rd = (sqrt(x_dest*x_dest + y_dest*y_dest)) / ((double*)params)[4]; // Normalized

	rs = rd;
	f = (((((double*)params)[3] * rs + ((double*)params)[2]) * rs +
		((double*)params)[1]) * rs + ((double*)params)[0]) * rs;

	while (abs(f - rd) > R_EPS && iter++ < MAXITER)
	{
		rs = rs - (f - rd) / (((4 * ((double*)params)[3] * rs + 3 * ((double*)params)[2]) * rs +
			2 * ((double*)params)[1]) * rs + ((double*)params)[0]);

		f = (((((double*)params)[3] * rs + ((double*)params)[2]) * rs +
			((double*)params)[1]) * rs + ((double*)params)[0]) * rs;
	}

	scale = (rd != 0.0) ? rs / rd : 1.0f;
	//      printf("scale = %lg iter = %d\n", scale,iter);

	*x_src = x_dest * scale;
	*y_src = y_dest * scale;
	return 1;
}

int inv_vertical(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double coefficients[4]

	register double rs, rd, f, scale;
	int iter = 0;

	rd = abs(y_dest) / ((double*)params)[4]; // Normalized

	rs = rd;
	f = (((((double*)params)[3] * rs + ((double*)params)[2]) * rs +
		((double*)params)[1]) * rs + ((double*)params)[0]) * rs;

	while (abs(f - rd) > R_EPS && iter++ < MAXITER)
	{
		rs = rs - (f - rd) / (((4 * ((double*)params)[3] * rs + 3 * ((double*)params)[2]) * rs +
			2 * ((double*)params)[1]) * rs + ((double*)params)[0]);

		f = (((((double*)params)[3] * rs + ((double*)params)[2]) * rs +
			((double*)params)[1]) * rs + ((double*)params)[0]) * rs;
	}

	scale = rs / rd;
	//      printf("scale = %lg iter = %d\n", scale,iter);

	*x_src = x_dest;
	*y_src = y_dest * scale;
	return 1;
}

int resize(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double scale_horizontal, double scale_vertical;

	*x_src = x_dest * var0;
	*y_src = y_dest * var1;
	return 1;
}

int shear(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double shear_horizontal, double shear_vertical;
	//	printf( "Entered shear function \n");


	*x_src = x_dest + var0 * y_dest;
	*y_src = y_dest + var1 * x_dest;
	return 1;
}

int shearInv(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double shear_horizontal, double shear_vertical;
	//	printf( "Entered shear inv function \n");


	*y_src = (y_dest - var1 * x_dest) / (1 - var1 * var0);
	*x_src = (x_dest - var0 * *y_src);
	return 1;
}

int horiz(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double horizontal shift

	*x_src = x_dest + shift;
	*y_src = y_dest;
	return 1;
}

int vert(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double vertical shift

	*x_src = x_dest;
	*y_src = y_dest + shift;
	return 1;
}


int radial(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double coefficients[4], scale, correction_radius

	register double r, scale;

	r = (sqrt(x_dest*x_dest + y_dest*y_dest)) / ((double*)params)[4];
	if (r < ((double*)params)[5])
	{
		scale = ((((double*)params)[3] * r + ((double*)params)[2]) * r +
			((double*)params)[1]) * r + ((double*)params)[0];
	}
	else
		scale = 1000.0;

	*x_src = x_dest * scale;
	*y_src = y_dest * scale;
	return 1;
}

int vertical(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double coefficients[4]

	register double r, scale;

	r = y_dest / ((double*)params)[4];

	if (r < 0.0) r = -r;

	scale = ((((double*)params)[3] * r + ((double*)params)[2]) * r +
		((double*)params)[1]) * r + ((double*)params)[0];

	*x_src = x_dest;
	*y_src = y_dest * scale;
	return 1;
}

int deregister(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params: double coefficients[4]

	register double r, scale;

	r = y_dest / ((double*)params)[4];

	if (r < 0.0) r = -r;

	scale = (((double*)params)[3] * r + ((double*)params)[2]) * r +
		((double*)params)[1];

	*x_src = x_dest + abs(y_dest) * scale;
	*y_src = y_dest;
	return 1;
}




int persp_sphere(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params :  double Matrix[3][3], double distanceparam

	register double theta, s, r;
	double v[3];

#if 0	// old
	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / *((double*)((void**)params)[1]);
	phi = atan2(y_dest, x_dest);

	v[0] = *((double*)((void**)params)[1]) * sin(theta) * cos(phi);
	v[1] = *((double*)((void**)params)[1]) * sin(theta) * sin(phi);
	v[2] = *((double*)((void**)params)[1]) * cos(theta);

	matrix_inv_mult((double(*)[3]) ((void**)params)[0], v);

	theta = atan2(sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
	phi = atan2(v[1], v[0]);
	*x_src = *((double*)((void**)params)[1]) * theta * cos(phi);
	*y_src = *((double*)((void**)params)[1]) * theta * sin(phi);
#endif

	r = sqrt(x_dest * x_dest + y_dest * y_dest);
	theta = r / *((double*)((void**)params)[1]);
	if (r == 0.0)
		s = 0.0;
	else
		s = sin(theta) / r;

	v[0] = s * x_dest;
	v[1] = s * y_dest;
	v[2] = cos(theta);

	matrix_inv_mult((double(*)[3]) ((void**)params)[0], v);

	r = sqrt(v[0] * v[0] + v[1] * v[1]);
	if (r == 0.0)
		theta = 0.0;
	else
		theta = *((double*)((void**)params)[1]) * atan2(r, v[2]) / r;
	*x_src = theta * v[0];
	*y_src = theta * v[1];

	return 1;
}

int persp_rect(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
	// params :  double Matrix[3][3], double distanceparam, double x-offset, double y-offset

	double v[3];

	v[0] = x_dest + *((double*)((void**)params)[2]);
	v[1] = y_dest + *((double*)((void**)params)[3]);
	v[2] = *((double*)((void**)params)[1]);

	matrix_inv_mult((double(*)[3]) ((void**)params)[0], v);

	*x_src = v[0] * *((double*)((void**)params)[1]) / v[2];
	*y_src = v[1] * *((double*)((void**)params)[1]) / v[2];
	return 1;
}



int rect_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{

	*x_src = distanceparam * tan(x_dest / distanceparam);
	*y_src = y_dest / cos(x_dest / distanceparam);
	return 1;
}

int pano_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	*x_src = distanceparam * atan(x_dest / distanceparam);
	*y_src = y_dest * cos(*x_src / distanceparam);
	return 1;
}

int rect_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double  phi, theta;

	phi = x_dest / distanceparam;
	theta = -y_dest / distanceparam + PI / 2.0;
	if (theta < 0)
	{
		theta = -theta;
		phi += PI;
	}
	if (theta > PI)
	{
		theta = PI - (theta - PI);
		phi += PI;
	}

#if 0
	v[2] = *((double*)params) * sin(theta) * cos(phi);   //  x' -> z
	v[0] = *((double*)params) * sin(theta) * sin(phi);	//  y' -> x
	v[1] = *((double*)params) * cos(theta);				//  z' -> y

	phi = atan2(v[1], v[0]);
	//  old:
	//	theta = atan2( sqrt( v[0]*v[0] + v[1]*v[1] ) , v[2] );
	//	rho = *((double*)params) * tan( theta );
	//  new:
	rho = *((double*)params) * sqrt(v[0] * v[0] + v[1] * v[1]) / v[2];
	*x_src = rho * cos(phi);
	*y_src = rho * sin(phi);
#endif
#if 1
	* x_src = distanceparam * tan(phi);
	*y_src = distanceparam / (tan(theta) * cos(phi));
#endif
	// normalize phi to be in the -PI, PI range
	while (phi <= -PI)
		phi += 2 * PI;
	while (phi > PI)
		phi -= 2 * PI;

	// check if the point is "in front" of the camera
	if (phi < -PI / 2.0 || phi > PI / 2.0) {
		// behind, transform considered invalid
		return 0;
	}
	else
		return 1;
	// normalize phi to be in the -PI, PI range
	while (phi <= -PI)
		phi += 2 * PI;
	while (phi > PI)
		phi -= 2 * PI;

	// check if the point is "in front" of the camera
	if (phi < -PI / 2.0 || phi > PI / 2.0) {
		// behind, transform considered invalid
		return 0;
	}
	else
		return 1;

}
// This is the cylindrical projection
int pano_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	*x_src = x_dest;
	*y_src = distanceparam * tan(y_dest / distanceparam);
	return 1;
}

int erect_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	*x_src = x_dest;
	*y_src = distanceparam * atan(y_dest / distanceparam);
	return 1;
}

/** convert from erect to lambert azimuthal */
int lambertazimuthal_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	double phi, lambda, k1;
	lambda = x_dest / distanceparam;
	phi = y_dest / distanceparam;

	if (abs(cos(phi) * cos(lambda) + 1.0) <= EPSLN) {
		*x_src = distanceparam * 2;
		*y_src = 0;
		return 0;
	}

	k1 = sqrt(2.0 / (1 + cos(phi) * cos(lambda)));

	*x_src = distanceparam * k1 * cos(phi) * sin(lambda);
	*y_src = distanceparam * k1 * sin(phi);

	return 1;
}

/** convert from lambert azimuthal to erect */
int erect_lambertazimuthal(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{

	double x, y, ro, c;

	x = x_dest / distanceparam;
	y = y_dest / distanceparam;

	assert(!isnan(x));
	assert(!isnan(y));

	if (fabs(x) > PI || fabs(y) > PI) {
		*y_src = 0;
		*x_src = 0;
		return 0;
	}

	ro = hypot(x, y);

	if (fabs(ro) <= EPSLN)
	{
		*y_src = 0;
		*x_src = 0;
		return 1;
	}

	c = 2 * asin(ro / 2.0);

	*y_src = distanceparam * asin((y * sin(c)) / ro);


	if (fabs(ro * cos(c)) <= EPSLN) {
		*x_src = 0;
		return 1;
	}

	*x_src = distanceparam * atan2(x * sin(c), (ro * cos(c)));

	return 1;
}

/** convert from erect to hammer */
int hammer_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	if (lambertazimuthal_erect(x_dest / 2.0, y_dest, x_src, y_src, params))
	{
		*x_src *= 2.0;
		return 1;
	}
	else
	{
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
}

/** convert from hammer to erect */
int erect_hammer(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	double x, y, z;
	x = x_dest / distanceparam;
	y = y_dest / distanceparam;
	z = 1.0 - (x * x / 16.0) - (y * y / 4.0);
	if (z<0)
	{
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	z = sqrt(z);
	*x_src = 2.0 * atan2(z*x, 2.0*(2.0*z*z - 1.0));
	*y_src = asin(y*z);
	if (fabs(*x_src) > PI || fabs(*y_src) > HALF_PI)
	{
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	*x_src *= distanceparam;
	*y_src *= distanceparam;
	return 1;
}

/** convert from erect to mercator FORWARD */
int mercator_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	*x_src = x_dest;
	*y_src = distanceparam*log(tan(y_dest / distanceparam) + 1 / cos(y_dest / distanceparam));
	return 1;
}

/** convert from mercator to erect INVERSE */
int erect_mercator(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	*x_src = x_dest;
	*y_src = distanceparam*atan(sinh(y_dest / distanceparam));
	return 1;
}


/** convert from erect to miller */
int millercylindrical_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	double phi, tanPhi;

	*x_src = x_dest;
	phi = y_dest / distanceparam;
	tanPhi = tan(PI / 4 + 0.4 * phi);
	if (tanPhi < 0) {
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	*y_src = distanceparam*log(tanPhi) / 0.8;
	return 1;
}

/** convert from erect to miller */
int arch_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	if (y_dest < 0) {
		return millercylindrical_erect(x_dest, y_dest, x_src, y_src, params);
	}
	else {
		return lambert_erect(x_dest, y_dest, x_src, y_src, params);
	}
}

/** convert from erect to miller */
int erect_arch(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	if (y_dest < 0) {
		return erect_millercylindrical(x_dest, y_dest, x_src, y_src, params);
	}
	else {
		return erect_lambert(x_dest, y_dest, x_src, y_src, params);
	}
}



/** convert from miller cylindrical to erect */
int erect_millercylindrical(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	double y;

	*x_src = x_dest;
	y = y_dest / distanceparam;
	y = 1.25 * atan(sinh(4 * y / 5.0));
	if (fabs(y) > HALF_PI) {
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	*y_src = distanceparam * y;
	return 1;
}


/** convert from erect to panini */
int panini_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	// this is the inverse

	double phi, lambdaHalf, temp, y, x;

	phi = y_dest / distanceparam;
	lambdaHalf = x_dest / (distanceparam * 2);
	x = 2 * tan(lambdaHalf);

	// Conver from central cylindrical
	phi = tan(phi);

	*x_src = distanceparam * x;
	temp = cos(lambdaHalf);

	y = tan(phi) / (temp * temp);

	// At this point we have mapped to central cylindrical, now to equirectangular

	*y_src = distanceparam *  y;

	return 1;
}

/** General Pannini Projection

setup_panini_general(&MakeParams) selects the Image struct
corresponding to the pannini_general image and returns its
address, or a NULL pointer for failure.

If the selected Image has an invalid precomputedCount, it
posts the distanceparam corresponding to min( max feasible
hFOV, requested hFOV) and puts working parameter values in
precomputeValue[] in the selected Image.

SetMakeParams (adjust.c) calls this function in lieu of setting
distanceparam.

The user-visible projection params, described in queryfeature.c,
are scaled to accomodate integer-valued control sliders in a GUI.
unscaleParams_panini_general() sets working values as follows:
cmpr 0:100:150 <-> d = 0:1:->infinity NOTE very nonlinear
tops, bots -100:100 <-> sqz -1:1 linear
< 0 gives soft squeeze
> 0 give transverse straightening squeeze
CAUTION these ranges are assumed, not read from queryfeature.c

maxFOVs_panini_general() calculates the maximum feasible FOVs
for a given scaled parameter set.   Those also depends on a
projection angle limit, that is hard coded here.  FOVs in degrees.

**/
#define MAX_PROJ_ANGLE 80

int unscaleParams_panini_general(
	double * gui_params,	// cmpr, tops, bots
	double * wrk_params		// d, t, b
	)
{
	double t;

	/* check for legal values */
	if (gui_params[0] < 0
		|| gui_params[0] > 150
		) return 0;
	if (gui_params[1] < -100
		|| gui_params[1] > 100
		) return 0;
	if (gui_params[2] < -100
		|| gui_params[2] > 100
		) return 0;

	/* post working param values */
	t = (150 - gui_params[0]) / 50;	/* 0:150 => 3:0 */
	wrk_params[0] = 1.5 / (t + 0.0001) - 1.5 / 3.0001;
	wrk_params[1] = gui_params[1] / 100;
	wrk_params[2] = gui_params[2] / 100;

	return 1;
}

/** convert from erect to lambert */
int lambert_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	*x_src = x_dest;
	*y_src = distanceparam*sin(y_dest / distanceparam);
	return 1;
}

/** convert from lambert to erect */
int erect_lambert(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	double y;
	*x_src = x_dest;
	y = y_dest / distanceparam;
	if (fabs(y) > 1) {
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	*y_src = distanceparam*asin(y);
	return 1;
}


/** convert from erect to transverse mercator */
int transmercator_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	// params: distanceparam
	double B;
	x_dest /= distanceparam;
	y_dest /= distanceparam;
	B = cos(y_dest)*sin(x_dest);
	*x_src = distanceparam * atanh(B);
	*y_src = distanceparam * atan2(tan(y_dest), cos(x_dest));

	if (isinf(*x_src)) {
		return 0;
	}

	return 1;
}

/** convert from erect to transverse mercator */
int erect_transmercator(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: distanceparam
	x_dest /= distanceparam;
	y_dest /= distanceparam;

	if (fabs(y_dest) > PI) {
		*y_src = 0;
		*x_src = 0;
		return 0;
	}


	*x_src = distanceparam * atan2(sinh(x_dest), cos(y_dest));
	*y_src = distanceparam * asin(sin(y_dest) / cosh(x_dest));
	return 1;
}

/** convert from erect to sinusoidal */
int sinusoidal_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	// params: distanceparam

	*x_src = distanceparam * (x_dest / distanceparam*cos(y_dest / distanceparam));
	*y_src = y_dest;
	return 1;
}

/** convert from sinusoidal to erect */
int erect_sinusoidal(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	// params: distanceparam

	*y_src = y_dest;
	*x_src = x_dest / cos(y_dest / distanceparam);
	if (*x_src / distanceparam < -PI || *x_src / distanceparam > PI)
		return 0;
	return 1;
}

/** convert from erect to stereographic */
int stereographic_erect_old(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	// params: distanceparam
	double lon = x_dest / distanceparam;
	double lat = y_dest / distanceparam;

	// use: R = 1
	double k = 2.0 / (1 + cos(lat)*cos(lon));
	*x_src = distanceparam * k*cos(lat)*sin(lon);
	*y_src = distanceparam * k*sin(lat);
	return 1;
}

int stereographic_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	double lon, lat;
	double sinphi, cosphi, coslon;
	double g, ksp;

	lon = x_dest / distanceparam;
	lat = y_dest / distanceparam;

	sinphi = sin(lat);
	cosphi = cos(lat);
	coslon = cos(lon);

	g = cosphi * coslon;

	// point projects to infinity:
	//    if (fabs(g + 1.0) <= EPSLN)

	ksp = distanceparam * 2.0 / (1.0 + g);
	*x_src = ksp * cosphi * sin(lon);
	*y_src = ksp * sinphi;

	return 1;
}

/** convert from stereographic to erect */
int erect_stereographic(double x_dest, double  y_dest, double* lon, double* lat, void*  params)
{
	double rh;		/* height above sphere*/
	double c;		/* angle					*/
	double sinc, cosc;	/* sin of c and cos of c			*/

	/* Inverse equations
	-----------------*/
	double x = x_dest / distanceparam;
	double y = y_dest / distanceparam;
	rh = sqrt(x * x + y * y);
	c = 2.0 * atan(rh / (2.0 * 1));
	sinc = sin(c);
	cosc = cos(c);
	*lon = 0;
	if (fabs(rh) <= EPSLN)
	{
		*lat = 0;
		return 0;
	}
	else
	{
		*lat = asin((y * sinc) / rh) * distanceparam;

		if ((fabs(cosc) < EPSLN) && (fabs(x) < EPSLN))
			return 0;
		else
			*lon = atan2((x * sinc), (cosc * rh)) * distanceparam;
	}
	return 1;
}


/** convert from stereographic to erect */
int erect_stereographic_old(double x_dest, double  y_dest, double* x_src, double* y_src, void*  params)
{
	// params: distanceparam

	// use: R = 1
	double p = sqrt(x_dest*x_dest + y_dest*y_dest) / distanceparam;
	double c = 2.0*atan(p / 2.0);

	*x_src = distanceparam * atan2(x_dest / distanceparam*sin(c), (p*cos(c)));
	*y_src = distanceparam * asin(y_dest / distanceparam*sin(c) / p);
	return 1;
}


int sphere_cp_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b

	register double phi, theta;
	phi = -x_dest / (var0 * PI / 2.0);
	theta = -(y_dest + var1) / (PI / 2.0);

	*x_src = theta * cos(phi);
	*y_src = theta * sin(phi);
	return 1;
}


int sphere_tp_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double phi, theta, r, s;
	double v[3];

	phi = x_dest / distanceparam;
	theta = -y_dest / distanceparam + PI / 2;
	if (theta < 0)
	{
		theta = -theta;
		phi += PI;
	}
	if (theta > PI)
	{
		theta = PI - (theta - PI);
		phi += PI;
	}


	s = sin(theta);
	v[0] = s * sin(phi);	//  y' -> x
	v[1] = cos(theta);				//  z' -> y

	r = sqrt(v[1] * v[1] + v[0] * v[0]);

	theta = distanceparam * atan2(r, s * cos(phi));

	*x_src = theta * v[0] / r;
	*y_src = theta * v[1] / r;
	return 1;
}


int erect_sphere_cp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b

	register double phi, theta;

#if 0
	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / var0;
	phi = atan2(y_dest, -x_dest);

	*x_src = var0 * phi;
	*y_src = var0 * theta - var1;
#endif
	theta = sqrt(x_dest * x_dest + y_dest * y_dest);
	phi = atan2(y_dest, -x_dest);

	*x_src = var0 * phi;
	*y_src = theta - var1;
	return 1;
}


int rect_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double rho, theta, r;

#if 0
	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / distanceparam;
	phi = atan2(y_dest, x_dest);

	if (theta > PI / 2.0 || theta < -PI / 2.0)
		theta = PI / 2.0;

	rho = distanceparam * tan(theta);

	*x_src = rho * cos(phi);
	*y_src = rho * sin(phi);
#endif
	r = sqrt(x_dest * x_dest + y_dest * y_dest);
	theta = r / distanceparam;

	if (theta >= PI / 2.0)
		rho = 1.6e16;
	else if (theta == 0.0)
		rho = 1.0;
	else
		rho = tan(theta) / theta;
	*x_src = rho * x_dest;
	*y_src = rho * y_dest;
	return 1;
}


int sphere_tp_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double  theta, r;

#if 0
	theta = atan(sqrt(x_dest*x_dest + y_dest*y_dest) / *((double*)params));
	phi = atan2(y_dest, x_dest);

	*x_src = *((double*)params) * theta * cos(phi);
	*y_src = *((double*)params) * theta * sin(phi);
#endif
	r = sqrt(x_dest*x_dest + y_dest*y_dest) / distanceparam;
	if (r == 0.0)
		theta = 1.0;
	else
		theta = atan(r) / r;

	*x_src = theta * x_dest;
	*y_src = theta * y_dest;
	return 1;
}


int sphere_tp_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double r, s, Phi, theta;

#if 0
	register double Theta, phi;
	double v[3];

	Phi = x_dest / *((double*)params);
	Theta = PI / 2.0 - atan(y_dest / distanceparam);


	v[2] = *((double*)params) * sin(Theta) * cos(Phi);   //  x' -> z
	v[0] = *((double*)params) * sin(Theta) * sin(Phi);	//  y' -> x
	v[1] = *((double*)params) * cos(Theta);				//  z' -> y

	theta = atan2(sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
	phi = atan2(v[1], v[0]);

	*x_src = *((double*)params) * theta * cos(phi);
	*y_src = *((double*)params) * theta * sin(phi);
#endif
#if 1
	Phi = x_dest / distanceparam;

	s = distanceparam * sin(Phi);	//  y' -> x

	r = sqrt(s*s + y_dest*y_dest);
	theta = distanceparam * atan2(r, (distanceparam * cos(Phi))) / r;

	*x_src = theta * s;
	*y_src = theta * y_dest;
#endif
	return 1;
}


int pano_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam
	register double r, s, theta;
	double v[3];

#if 0
	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / distanceparam;
	phi = atan2(y_dest, x_dest);

	v[1] = *((double*)params) * sin(theta) * cos(phi);   //  x' -> y
	v[2] = *((double*)params) * sin(theta) * sin(phi);	//  y' -> z
	v[0] = *((double*)params) * cos(theta);				//  z' -> x

	theta = atan2(sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
	phi = atan2(v[1], v[0]);

	*x_src = *((double*)params) * phi;
	*y_src = *((double*)params) * tan((-theta + PI / 2.0));
#endif

	r = sqrt(x_dest * x_dest + y_dest * y_dest);
	theta = r / distanceparam;
	if (theta == 0.0)
		s = 1.0 / distanceparam;
	else
		s = sin(theta) / r;

	v[1] = s * x_dest;   //  x' -> y
	v[0] = cos(theta);				//  z' -> x


	*x_src = distanceparam * atan2(v[1], v[0]);
	*y_src = distanceparam * s * y_dest / sqrt(v[0] * v[0] + v[1] * v[1]);

	return 1;
}


int sphere_cp_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double phi, theta;


	phi = -x_dest / (distanceparam * PI / 2.0);
	theta = PI / 2.0 + atan(y_dest / (distanceparam * PI / 2.0));

	*x_src = distanceparam * theta * cos(phi);
	*y_src = distanceparam * theta * sin(phi);
	return 1;
}


int erect_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam
#if 0
	theta = atan(sqrt(x_dest*x_dest + y_dest*y_dest) / distanceparam);
	phi = atan2(y_dest, x_dest);


	v[1] = distanceparam * sin(theta) * cos(phi);   //  x' -> y
	v[2] = distanceparam * sin(theta) * sin(phi);	//  y' -> z
	v[0] = distanceparam * cos(theta);				//  z' -> x

	theta = atan2(sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
	phi = atan2(v[1], v[0]);

	*x_src = distanceparam * phi;
	*y_src = distanceparam * (-theta + PI / 2.0);
#endif

	*x_src = distanceparam * atan2(x_dest, distanceparam);
	*y_src = distanceparam * atan2(y_dest, sqrt(distanceparam*distanceparam + x_dest*x_dest));

	return 1;
}

/** convert erect to cartesian XYZ coordinates
*/
int cart_erect(double x_dest, double y_dest, double * xyz, double distance)
{
	// phi is azimuth (negative angle around y axis, starting at the z axis)
	double phi = x_dest / distance;
	double theta_zenith = PI / 2.0 - (y_dest / distance);
	// compute cartesian coordinates..
	//pos[2] = cos(-phi)*sin(theta_zenith);
	//pos[0] = sin(-phi)*sin(theta_zenith);
	//pos[1] = cos(theta_zenith);

	xyz[0] = sin(theta_zenith)*sin(phi);
	xyz[1] = cos(theta_zenith);
	xyz[2] = sin(theta_zenith)*-cos(phi);

	return 1;
}


/** convert cartesian coordinates into spherical ones
*/
int erect_cart(double * xyz, double *x_src, double *y_src, double distance)
{
	*x_src = atan2(xyz[0], -xyz[2]) * distance;
	*y_src = asin(xyz[1] / sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2])) * distance;

	return 1;
}


/** Compute intersection between line and point.
*  n : a,b,c,d coefficients of plane (a,b,c = normal vector)
*  p1: point on line
*  p2: point on line
*  See http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
*/
int line_plane_intersection(double n[4],
	double p1[3],
	double p2[3],
	double * result)
{
	int i;
	// direction vector of line
	double d[3];
	double u, num, den;

	for (i = 0; i<3; i++)
		d[i] = p2[i] - p1[i];
	num = n[0] * p1[0] + n[1] * p1[1] + n[2] * p1[2] + n[3];
	den = -n[0] * d[0] - n[1] * d[1] - n[2] * d[2];
	if (fabs(den) < 1e-15) {
		return 0;
	}
	u = num / den;

	if (u < 0) {
		// This is match is in the wrong direction, ignore
		return 0;
	}
	/* printf("intersect, dir: %f %f %f, num: %f, denom: %f, u: %f\n", d[0], d[1], d[2], num, den, u);
	*/

	for (i = 0; i<3; i++)
		result[i] = p1[i] + u*d[i];

	return 1;
}

int erect_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam

	register double  theta, r, s;
	double	v[3];
#if 0
	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / *((double*)params);
	phi = atan2(y_dest, x_dest);

	v[1] = *((double*)params) * sin(theta) * cos(phi);   //  x' -> y
	v[2] = *((double*)params) * sin(theta) * sin(phi);	//  y' -> z
	v[0] = *((double*)params) * cos(theta);				//  z' -> x

	theta = atan(sqrt(v[0] * v[0] + v[1] * v[1]) / v[2]); //was atan2
	phi = atan2(v[1], v[0]);

	*x_src = *((double*)params) * phi;
	if (theta > 0.0)
	{
		*y_src = *((double*)params) * (-theta + PI / 2.0);
	}
	else
		*y_src = *((double*)params) * (-theta - PI / 2.0);
#endif

	r = sqrt(x_dest * x_dest + y_dest * y_dest);
	theta = r / distanceparam;
	if (theta == 0.0)
		s = 1.0 / distanceparam;
	else
		s = sin(theta) / r;

	v[1] = s * x_dest;
	v[0] = cos(theta);


	*x_src = distanceparam * atan2(v[1], v[0]);
	*y_src = distanceparam * atan(s * y_dest / sqrt(v[0] * v[0] + v[1] * v[1]));
	return 1;
}


int mirror_sphere_cp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b

	register double rho, phi, theta;

	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / ((double*)params)[0];
	phi = atan2(y_dest, x_dest);

	rho = ((double*)params)[1] * sin(theta / 2.0);

	*x_src = -rho * cos(phi);
	*y_src = rho * sin(phi);
	return 1;
}


int mirror_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b, double b2

	register double phi, theta, rho;

	phi = x_dest / (((double*)params)[0] * PI / 2.0);
	theta = -(y_dest + ((double*)params)[2]) / (((double*)params)[0] * PI / 2.0);

	rho = ((double*)params)[1] * sin(theta / 2.0);

	*x_src = -rho * cos(phi);
	*y_src = rho * sin(phi);
	return 1;
}


int mirror_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b

	register double phi, theta, rho;


	phi = -x_dest / (((double*)params)[0] * PI / 2.0);
	theta = PI / 2.0 + atan(y_dest / (((double*)params)[0] * PI / 2.0));

	rho = ((double*)params)[1] * sin(theta / 2.0);

	*x_src = rho * cos(phi);
	*y_src = rho * sin(phi);
	return 1;
}


int sphere_cp_mirror(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distanceparam, double b

	register double phi, theta, rho;

	rho = sqrt(x_dest*x_dest + y_dest*y_dest);

	theta = 2 * asin(rho / ((double*)params)[1]);
	phi = atan2(y_dest, x_dest);

	*x_src = ((double*)params)[0] * theta * cos(phi);
	*y_src = ((double*)params)[0] * theta * sin(phi);
	return 1;
}

int sphere_tp_mirror(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	register double c;
	register double normalizedX, normalizedY;
	register double azi;

	normalizedX = x_dest / distanceparam;
	normalizedY = y_dest / distanceparam;

	c = 2 * asin(hypot(normalizedX, normalizedY));
	azi = atan2(y_dest, x_dest);

	*x_src = distanceparam * c * cos(azi);
	*y_src = distanceparam * c * sin(azi);

	return 1;
}

int mirror_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{

	register double c;
	register double normalizedX, normalizedY;
	register double azi;
	normalizedY = y_dest / distanceparam;
	normalizedX = x_dest / distanceparam;


	c = hypot(normalizedX, normalizedY);
	azi = atan2(y_dest, x_dest);

	*x_src = distanceparam * sin(c / 2) * cos(azi);
	*y_src = distanceparam * sin(c / 2) * sin(azi);

	return 1;
}

int sphere_tp_equisolid(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distance

	register double phi, theta, rho;

	rho = sqrt(x_dest*x_dest + y_dest*y_dest);

	theta = 2.0 * asin(rho / (2.0*((double*)params)[0]));
	phi = atan2(y_dest, x_dest);

	*x_src = ((double*)params)[0] * theta * cos(phi);
	*y_src = ((double*)params)[0] * theta * sin(phi);
	return 1;
}


int equisolid_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distance

	register double rho, phi, theta;

	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / ((double*)params)[0];
	phi = atan2(y_dest, x_dest);

	rho = 2.0 * ((double*)params)[0] * sin(theta / 2.0);

	*x_src = rho * cos(phi);
	*y_src = rho * sin(phi);
	return 1;
}

int sphere_tp_orthographic(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distance

	register double phi, theta, rho;

	rho = sqrt(x_dest*x_dest + y_dest*y_dest);

	// orthographic projection is limited to fov of 180 deg
	if (rho>((double*)params)[0])
	{
		*x_src = 0;
		*y_src = 0;
		return 0;
	};
	theta = 1.0 * asin(rho / (1.0*((double*)params)[0]));
	phi = atan2(y_dest, x_dest);


	*x_src = ((double*)params)[0] * theta * cos(phi);
	*y_src = ((double*)params)[0] * theta * sin(phi);

	return 1;
}

int orthographic_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params)
{
	// params: double distance

	register double rho, phi, theta;

	theta = sqrt(x_dest * x_dest + y_dest * y_dest) / ((double*)params)[0];
	phi = atan2(y_dest, x_dest);
	//orthographic projection is limited to fov of 180 deg
	if (fabs(theta)>HALF_PI)
	{
		*x_src = 0;
		*y_src = 0;
		return 0;
	};

	rho = 1.0 * ((double*)params)[0] * sin(theta / 1.0);

	*x_src = rho * cos(phi);
	*y_src = rho * sin(phi);
	return 1;
}

int shift_scale_rotate(double x_dest, double  y_dest, double* x_src, double* y_src, void* params){
	// params: double shift_x, shift_y, scale, cos_phi, sin_phi

	register double x = x_dest - ((double*)params)[0];
	register double y = y_dest - ((double*)params)[1];

	*x_src = (x * ((double*)params)[3] - y * ((double*)params)[4]) * ((double*)params)[2];
	*y_src = (x * ((double*)params)[4] + y * ((double*)params)[3]) * ((double*)params)[2];
	return 1;
}


void CubeZero(double *a, int *n, double *root) {
	if (a[3] == 0.0){ // second order polynomial
		SquareZero(a, n, root);
	}
	else{
		double p = ((-1.0 / 3.0) * (a[2] / a[3]) * (a[2] / a[3]) + a[1] / a[3]) / 3.0;
		double q = ((2.0 / 27.0) * (a[2] / a[3]) * (a[2] / a[3]) * (a[2] / a[3]) - (1.0 / 3.0) * (a[2] / a[3]) * (a[1] / a[3]) + a[0] / a[3]) / 2.0;

		if (q*q + p*p*p >= 0.0){
			*n = 1;
			root[0] = CubeRoot(-q + sqrt(q*q + p*p*p)) + CubeRoot(-q - sqrt(q*q + p*p*p)) - a[2] / (3.0 * a[3]);
		}
		else{
			double phi = acos(-q / sqrt(-p*p*p));
			*n = 3;
			root[0] = 2.0 * sqrt(-p) * cos(phi / 3.0) - a[2] / (3.0 * a[3]);
			root[1] = -2.0 * sqrt(-p) * cos(phi / 3.0 + PI / 3.0) - a[2] / (3.0 * a[3]);
			root[2] = -2.0 * sqrt(-p) * cos(phi / 3.0 - PI / 3.0) - a[2] / (3.0 * a[3]);
		}
	}
	// PrintError("%lg, %lg, %lg, %lg root = %lg", a[3], a[2], a[1], a[0], root[0]);
}

void SquareZero(double *a, int *n, double *root){
	if (a[2] == 0.0){ // linear equation
		if (a[1] == 0.0){ // constant
			if (a[0] == 0.0){
				*n = 1; root[0] = 0.0;
			}
			else{
				*n = 0;
			}
		}
		else{
			*n = 1; root[0] = -a[0] / a[1];
		}
	}
	else{
		if (4.0 * a[2] * a[0] > a[1] * a[1]){
			*n = 0;
		}
		else{
			*n = 2;
			root[0] = (-a[1] + sqrt(a[1] * a[1] - 4.0 * a[2] * a[0])) / (2.0 * a[2]);
			root[1] = (-a[1] - sqrt(a[1] * a[1] - 4.0 * a[2] * a[0])) / (2.0 * a[2]);
		}
	}

}

double CubeRoot(double x){
	if (x == 0.0)
		return 0.0;
	else if (x > 0.0)
		return pow(x, 1.0 / 3.0);
	else
		return -pow(-x, 1.0 / 3.0);
}

double SmallestRoot(double *p){
	int n, i;
	double root[3], sroot = 1000.0;

	CubeZero(p, &n, root);

	for (i = 0; i<n; i++){
		// PrintError("Root %d = %lg", i,root[i]);
		if (root[i] > 0.0 && root[i] < sroot)
			sroot = root[i];
	}

	// PrintError("Smallest Root  = %lg", sroot);
	return sroot;
}

void SetCorrectionRadius(double* rad) {
	double a[4];
	int i, k;

	for (i = 0; i<3; i++)
	{
		for (k = 0; k<4; k++)
		{
			a[k] = 0.0;//1.0e-10;
			if (rad[k] != 0.0)
			{
				a[k] = (k + 1) * rad[k];
			}
		}
		rad[4] = SmallestRoot(a);
	}
}

}
