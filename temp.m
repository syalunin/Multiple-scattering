clear
close
clc
ry = 1/2;
ry2ev = 13.6056919;
ev = 1/27.2113838;
ev2kelv = 11605;
kelv = ev/ev2kelv;
ry2kelv = ry2ev*ev2kelv;
energy = 86*ev;
velocity = sqrt(2*energy);

% incident wave
ki = velocity*(1+0.1i);
theta = pi;
phi = 0;
[ki_dir(1),ki_dir(2),ki_dir(3)] = sph2cart(phi,pi/2-theta,1);
ki_dir = reshape(ki_dir,3,1);
nr = 3;

load('parameters_100_50.mat');
load('Gaunt_factors.mat');
fprintf(' parameters loaded \n');

% Compute matrix elements by summing the scattering contributions
% from neighbouring Wigner–Seitz cells, 2*nr+1 in every direction

% loops over adjacent Wigner–Seitz cells
tic
H = zeros(nat,nat,(2*lmax+1)^2);
for l=0:2*lmax
  for m=-l:l
    n_lm = l*(l+1)+m+1;
    for na=1:nat
      for nb=1:nat
        for n1=-nr:nr
          for n2=-nr:nr
            for n3=0
              if n1^2+n2^2+n3^2~=0
                r = n1*at(:,1)+n2*at(:,2)+n3*at(:,3);
                exp_kr = exp(1i*ki*ki_dir'*r);
                r = r+tau(:,na)-tau(:,nb);
                [phi,elev,r] = cart2sph(r(1),r(2),r(3));
                H(na,nb,n_lm) = 1i^l*sqrt(pi/2/ki/r)*besselh(l+1/2,ki*r)*...
                                spherical_harmonic(l,m,pi/2-elev,phi)*exp_kr;
              end
            end
          end
        end
      end
    end
  end
end

% compute matrix g
A = zeros((lmax+1)^2*nat);
count1 = 0;
for na=1:nat
  for l=0:lmax
    for m=-l:l
      count1 = count1+1;
      count2 = 0;
      for nb=1:nat
        for ll=0:lmax
          for mm=-ll:ll
            count2 = count2+1;
            if ll<l
              l1 = l; l2 = ll;
              m1 = m; m2 = -mm;
            else
              l2 = l; l1 = ll;
              m2 = m; m1 = -mm;
            end
            if m2<0
              m1 = -m1;
              m2 = -m2;
            end
            n = l1*(l1+1)*(l1+2)*(3*l1-1)/12+(l1+1)*(l1+2)*(l1+m1)/2+l2*(l2+1)/2+m2+1;
            count = index(n);
            m3 = -m1-m2;
            l3_min = max(l1-l2,abs(m3));
            for lll=l1+l2:-2:l3_min
              n_lm = lll*(lll+1)+mm-m+1;
              A(count1,count2) = A(count1,count2)+4*pi*(-1)^mm*Gaunt(count)*H(na,nb,n_lm);
              count = count+1;
            end
          end
        end
      end
    end
  end
end
toc
