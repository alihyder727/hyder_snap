#include "../athena.hpp"

inline Real dfdx(AthenaArray<Real> const& f, Real dx, int i, int j, int k)
{
  return (f(k+1,j+1,i+1) + f(k+1,j,i+1) + f(k,j+1,i+1) + f(k,j,i+1) -
          f(k+1,j+1,i  ) - f(k+1,j,i  ) - f(k,j+1,i  ) - f(k,j,i  ))/(4*dx);
}

inline Real dfdy(AthenaArray<Real> const& f, Real dy, int i, int j, int k)
{
  return (f(k+1,j+1,i+1) + f(k+1,j+1,i) + f(k,j+1,i+1) + f(k,j+1,i) -
          f(k+1,j  ,i  ) - f(k+1,j  ,i) - f(k,j  ,i  ) - f(k,j  ,i))/(4*dy);
}

inline Real dfdz(AthenaArray<Real> const& f, Real dz, int i, int j, int k)
{
  return (f(k+1,j+1,i+1) + f(k+1,j,i+1) + f(k+1,j+1,i) + f(k+1,j,i) -
          f(k  ,j+1,i+1) - f(k  ,j,i+1) - f(k  ,j+1,i) - f(k  ,j,i))/(4*dz);
}

inline Real d2fdxy(AthenaArray<Real> const& f, Real dx, Real dy, int i, int j, int k)
{
  return (f(k  ,j+1,i+1) + f(k  ,j,i) - f(k  ,j+1,i) - f(k  ,j,i+1) +
          f(k+1,j+1,i+1) + f(k+1,j,i) - f(k+1,j+1,i) - f(k+1,j,i+1))/(2*dx*dy);
}

inline Real d2fdyz(AthenaArray<Real> const& f, Real dy, Real dz, int i, int j, int k)
{
  return (f(k+1,j+1,i  ) + f(k,j,i  ) - f(k+1,j,i  ) - f(k,j+1,i  ) +
          f(k+1,j+1,i+1) + f(k,j,i+1) - f(k+1,j,i+1) - f(k,j+1,i+1))/(2*dy*dz);
}

inline Real d2fdzx(AthenaArray<Real> const& f, Real dz, Real dx, int i, int j, int k)
{
  return (f(k+1,j  ,i+1) + f(k,j  ,i) - f(k+1,j  ,i) - f(k,j  ,i+1) +
          f(k+1,j+1,i+1) + f(k,j+1,i) - f(k+1,j+1,i) - f(k,j+1,i+1))/(2*dz*dx);
}

inline Real d2fdxx(AthenaArray<Real> const& f, Real dx, int i, int j, int k)
{
  return (f(k  ,j+1,i+2) - f(k  ,j+1,i+1) - f(k  ,j+1,i) + f(k  ,j+1,i-1) +
          f(k  ,j  ,i+2) - f(k  ,j  ,i+1) - f(k  ,j  ,i) + f(k  ,j,  i-1) +
          f(k+1,j+1,i+2) - f(k+1,j+1,i+1) - f(k+1,j+1,i) + f(k+1,j+1,i-1) +
          f(k+1,j  ,i+2) - f(k+1,j  ,i+1) - f(k+1,j  ,i) + f(k+1,j,  i-1))/(8*dx*dx);
}

inline Real d2fdyy(AthenaArray<Real> const& f, Real dy, int i, int j, int k)
{
  return (f(k  ,j+2,i+1) - f(k  ,j+1,i+1) - f(k  ,j,i+1) + f(k  ,j-1,i+1) +
          f(k  ,j+2,i  ) - f(k  ,j+1,i  ) - f(k  ,j,i  ) + f(k  ,j-1,i  ) + 
          f(k+1,j+2,i+1) - f(k+1,j+1,i+1) - f(k+1,j,i+1) + f(k+1,j-1,i+1) +
          f(k+1,j+2,i  ) - f(k+1,j+1,i  ) - f(k+1,j,i  ) + f(k+1,j-1,i  ))/(8*dy*dy);
}

inline Real d2fdzz(AthenaArray<Real> const& f, Real dz, int i, int j, int k)
{
  return (f(k+2,j  ,i+1) - f(k+1,j  ,i+1) - f(k,j  ,i+1) + f(k-1,j  ,i+1) +
          f(k+2,j  ,i  ) - f(k+1,j  ,i  ) - f(k,j  ,i  ) + f(k-1,j  ,i  ) +
          f(k+2,j+1,i+1) - f(k+1,j+1,i+1) - f(k,j+1,i+1) + f(k-1,j+1,i+1) +
          f(k+2,j+1,i  ) - f(k+1,j+1,i  ) - f(k,j+1,i  ) + f(k-1,j+1,i  ))/(8*dz*dz);
}
