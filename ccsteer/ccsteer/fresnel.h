#ifndef FRESNEL_H
#define FRESNEL_H
#ifdef __cplusplus
extern "C" {
#endif

	/* fresnel_c(x) - Fresnel Cosine Integral
	* C(x)=fresnel_c(x)=\dint\limits_{0}^{x}\cos (\frac{\pi}{2}t^{2})dt
	*/
	double fresnel_c(double x);
	/* fresnel_s(x) - Fresnel Sine Integral
	* S(x)=fresnel_s(x)=\dint\limits_{0}^{x}\sin (\frac{\pi}{2}t^{2})dt
	*/
	double fresnel_s(double x);

#ifdef __cplusplus
}
#endif
#endif