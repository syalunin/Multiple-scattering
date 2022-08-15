function [rad,dens,pot] = denspot(Z,E)
% reads electron density distribution from the database
% note that all variables (excepting E) are given in atomic units,
% E should be in eV

str = sprintf("%0.3i",Z);
f = fopen("database/Z_"+str+".den",'r');
data = fscanf(f,'%s',10);
data = fscanf(f,'%f');
dens = data(2:2:end);
rad = data(1:2:end-1);
energy = E/27.2113838;

Vst = 0*rad;
for j=1:numel(rad)-1
  Vst(j) = 4*pi*integral(@(r) spline(rad,dens,r).*(1-r/rad(j)).*r,rad(j),rad(end));
end
Vex = (energy-Vst)/2-sqrt((energy-Vst).^2+4*pi*dens)/2;
pot = Vst+Vex;
fclose(f);
end

