#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_
#include "../defs.hpp"

#define SQR(x) ( (x)*(x) )

/* 2-nd order plm for non-uniform grid
template<typename T>
inline T interp_plm(T phim1, T phi, T phip1, Real dxl, Real dx) {
  Real dwl = phi - phim1;
  Real dwr = phip1 - phi;
  Real dw2 = dwl*dwr;
  Real dwm = 2.0*dw2/(dwl + dwr);
  if (dw2 <= 0.0) dwm = 0.0;
  return phi - dxl/dx*dwm;
}*/

// limiter
template<typename T>
inline T minmod(T a, T b) {
  return a*std::max(0., std::min(1., b/a));
}

template<typename T>
inline T superbee(T a, T b) {
  T r = b/a;
  return a*std::max(0., std::min(1., 2.*r), std::min(2., r));
}

template<typename T>
inline T vanleer(T a, T b) {
  T r = b/a;
  return a*(r + fabs(r))/(1. + fabs(r));
}

template<typename T>
inline T mclimiter(T a, T b) {
  T r = b/a;
  T c = (1. + r)/2.;
  return std::max(0., std::min(std::min(c, 2.), 2.*r));
}

// sign
template <typename T> 
inline int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

// 2-rd polynomial 
template<typename T>
inline T interp_cp2(T phim1, T phi) {
  return 0.5*phim1 + 0.5*phi;
}

// 1-rd upwind-biased polynomial 
template<typename T>
inline T interp_bp2(T phim1, T phi, int sgn) {
  T phih = interp_cp2(phim1, phi);
  T phid = 0.5*(phi - phim1);
  return phih - sgn*phid;
}

#ifdef UNIFORM_GRID

// 3-th polynomial 
template<typename T>
inline T interp_cp3(T phim1, T phi, T phip1) {
  return 1./3.*phim1 + 5./6.*phi - 1./6.*phip1;
};

// 4-th polynomial 
template<typename T>
inline T interp_cp4(T phim2, T phim1, T phi, T phip1) {
  return -1./12.*phim2 + 7./12.*phim1 + 7./12.*phi - 1./12.*phip1;
};

// 3-rd upwind-biased polynomial 
template<typename T>
inline T interp_bp3(T phim2, T phim1, T phi, T phip1, int sgn) {
  T phih = interp_cp4(phim2, phim1, phi, phip1);
  T phid = ((phip1 - phim2) - 3.*(phi - phim1))/12.;
  return phih + sgn*phid;
}

// WENO 3 interpolation
template<typename T>
inline T interp_weno3(T phim1, T phi, T phip1) {
  T p0 = (1.0/2.0)*phi + (1.0/2.0)*phim1;
  T p1 = (-1.0/2.0)*phip1 + (3.0/2.0)*phi;

  T beta0 = (phim1 - phi)*(phim1 - phi);
  T beta1 = (phi - phip1)*(phi - phip1);

  T alpha0 = (2.0/3.0)/SQR(beta0 + 1e-10);
  T alpha1 = (1.0/3.0)/SQR(beta1 + 1e-10);

  return (alpha0*p0 + alpha1*p1)/(alpha0 + alpha1);
};

// 5-th polynomial 
template<typename T>
inline T interp_cp5(T phim2, T phim1, T phi, T phip1, T phip2) {
  return -1./20.*phim2 + 9./20.*phim1 + 47./60.*phi - 13./60.*phip1 + 1./30.*phip2;
};

// 6-th polynomial 
template<typename T>
inline T interp_cp6(T phim3, T phim2, T phim1, T phi, T phip1, T phip2) {
  return 1./60.*phim3 - 2./15.*phim2 + 37./60.*phim1 + 37./60.*phi - 2./15.*phip1 + 1./60.*phip2;
};

// 5-th upwind-biased polynomial 
template<typename T>
inline T interp_bp5(T phim3, T phim2, T phim1, T phi, T phip1, T phip2, int sgn) {
  T phih = interp_cp6(phim3, phim2, phim1, phi, phip1, phip2);
  T phid = ((phip2 - phim3) - 5.*(phip1 - phim2) + 10.*(phi - phim1))/60.;
  return phih - sgn*phid;
};

// WENO 5 interpolation
template<typename T>
inline T interp_weno5(T phim2, T phim1, T phi, T phip1, T phip2) {
  T p0 = (1./3.)*phi + (5./6.)*phim1 - (1./6.)*phim2;
  T p1 = (-1./6.)*phip1 + (5./6.)*phi + (1./3.)*phim1;
  T p2 = (1./3.)*phip2 - (7./6.)*phip1 + (11./6.)*phi;

  T beta0 = 13./12.*SQR(phi - 2.*phim1 + phim2) + .25*SQR(3.*phi - 4.*phim1 + phim2);
  T beta1 = 13./12.*SQR(phip1 - 2.*phi + phim1) + .25*SQR(phip1 - phim1);
  T beta2 = 13./12.*SQR(phip2 - 2.*phip1 + phi) + .25*SQR(phip2 - 4.*phip1 + 3.*phi);

  T alpha0 = .3/SQR(beta0 + 1e-10);
  T alpha1 = .6/SQR(beta1 + 1e-10);
  T alpha2 = .1/SQR(beta2 + 1e-10);

  return (alpha0*p0 + alpha1*p1 + alpha2*p2)/(alpha0 + alpha1 + alpha2);
};

#endif  // UNIFORM_GRID

// 3rd order polynomial with inflection point
template<typename T>
inline T inflection3_cell1(T f1, T f2, T f3) {
  return 7./3.*f1 - 5./3.*f2 + 1./3.*f3;
}

template<typename T>
inline T inflection3_cell2(T f1, T f2, T f3) {
  return 10./3.*f1 - 8./3.*f2 + 1./3.*f3;
}

template<typename T>
inline T inflection3_cell3(T f1, T f2, T f3) {
  return 10./3.*f1 - 5./3.*f2 - 2./3.*f3;
}

#undef SQR

#endif
