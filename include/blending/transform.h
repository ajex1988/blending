#ifndef INCLUDE_BLENDING_TRANSFORM_H_
#define INCLUDE_BLENDING_TRANSFORM_H_

namespace blending {

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef bzero
#define bzero(dest, len)   memset((dest), 0, (len))
#endif


//---------------------- Some useful math defines --------------------------

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif
#ifndef HALF_PI
#define HALF_PI (PI*0.5)
#endif

#define EPSLN	1.0e-10

// Normalize an angle to +/-180degrees

#define NORM_ANGLE( x )      while( x >180.0 ) x -= 360.0; while( x < -180.0 ) x += 360.0;
#define NORM_ANGLE_RAD( x )  while( (x) >PI ) (x) -= 2 * PI; while( (x) < -PI ) (x) += 2 * PI;

// Convert degree to radian

#define DEG_TO_RAD( x )		( (x) * 2.0 * PI / 360.0 )

// and reverse

#define RAD_TO_DEG( x )		( (x) * 360.0 / ( 2.0 * PI ) )

// Convert double x to unsigned char/short c



#define	DBL_TO_UC( c, x )	if((x)>255.0) c=255U;								\
								else if ((x)<0.0) c = 0;							\
								else c = (unsigned char)floor((x)+0.5);

#define	DBL_TO_US( c, x )	if((x)>65535.0) c=65535U;							\
								else if ((x)<0.0) c = 0;							\
								else c = (unsigned short)floor((x)+0.5);

#define	DBL_TO_FL( c, x )	if((x)>1e+038) c=1e+038;							\
								else if ((x)<0.0) c = 0;							\
								else c = (float)(x);

void 	SetMatrix(double a, double b, double c, double m[3][3], int cl);

int resize(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int shear(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int shearInv(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int horiz(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int vert(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int radial(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int radial_brown(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int persp_sphere(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int persp_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int rect_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int pano_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int pano_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_cp_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_tp_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_sphere_cp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int rect_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_tp_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_cp_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int rect_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_rect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int plane_transfer_to_camera(double x_dest, double y_dest, double * x_src, double * y_src, void * params);
int plane_transfer_from_camera(double x_dest, double y_dest, double * x_src, double * y_src, void * params);
int erect_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int mirror_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int mercator_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_mercator(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int lambert_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_lambert(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_lambertazimuthal(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int lambertazimuthal_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_hammer(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int hammer_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int transmercator_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_transmercator(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sinusoidal_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_sinusoidal(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int stereographic_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_stereographic(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int albersequalareaconic_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_albersequalareaconic(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int albersequalareaconic_distance(double *x_src, void* params);
int millercylindrical_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_millercylindrical(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int panini_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_panini(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int equipanini_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_equipanini(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int panini_general_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_panini_general(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int arch_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_arch(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int biplane_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_biplane(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int biplane_distance(double width, double b, void* params);
int triplane_erect(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int erect_triplane(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int triplane_distance(double width, double b, void* params);

int mirror_sphere_cp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int mirror_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_cp_mirror(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int sphere_tp_pano(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int pano_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_tp_mirror(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int mirror_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_tp_equisolid(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int equisolid_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int sphere_tp_orthographic(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int orthographic_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int sphere_tp_thoby(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);
int thoby_sphere_tp(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);


int rotate_erect(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int inv_radial(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int inv_radial_brown(double x_dest, double y_dest, double* x_src, double* y_src, void* params);

int vertical(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int inv_vertical(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int deregister(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int tmorph(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);

int shift_scale_rotate(double x_dest, double  y_dest, double* x_src, double* y_src, void* params);


void SetCorrectionRadius(double* rad);

}
#endif // INCLUDE_BLENDING_TRANSFORM_H_
