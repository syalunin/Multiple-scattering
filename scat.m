function [amp,del,theta] = scat(E,lmax,rad,pot)
% Computes scattering amplitudes and phase shifts with l<=lmax
% note that all variables (excepting E) are given in atomic units,
% E should be specified in eV

energy = E/27.2113838;
k = sqrt(2*energy);
z = k*rad(end);
fun = spline(rad,pot);
cf = fun.coefs;

del = zeros(lmax+1,1);
for l=0:lmax
  icond = [1 0];
  for n=1:numel(rad)-1
    rhs = @(r,sol) [2*(l*(l+1)/2./r.^2+cf(n,1)*(r-rad(n)).^3+...
    cf(n,2)*(r-rad(n)).^2+cf(n,3)*(r-rad(n))+cf(n,4)-energy)*sol(2) sol(1)]';
    [r,sol] = ode45(rhs,[rad(n) rad(n+1)],icond);
    icond = sol(end,:);
  end
  dpsi = icond(1);
  psi = icond(2);
  jl = sqrt(z)*besselj(l+1/2,z);
  yl = sqrt(z)*bessely(l+1/2,z);
  djl = sqrt(z)*(besselj(l-1/2,z)-l/z*besselj(l+1/2,z));
  dyl = sqrt(z)*(bessely(l-1/2,z)-l/z*bessely(l+1/2,z));
  del(l+1) = atan((dpsi*jl-k*psi*djl)/(dpsi*yl-k*psi*dyl));
end

amp = zeros(300,1);
theta = linspace(0,180,300)';
for l=0:lmax
  amp = amp+(2*l+1)*(exp(2i*del(l+1))-1)*legendreP(l,cos(pi*theta/180));
end
amp = amp/2i/k;
end
