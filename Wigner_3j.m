function result = Wigner_3j(j1,j2,j3,m1,m2,m3,fact)
% Compute 3j symbol using Racah formula and floating-point arithmetic

% 19-July-2022

result = 0;
if j3 < abs(j1-j2) || (j1+j2) < j3 || j1 < 0 || j2 < 0 || j3 < 0, return; end
if (m1+m2+m3) ~= 0 || j1 < abs(m1) || j2 < abs(m2) || j3 < abs(m3), return; end

kmin = max( 0, max(j2-j3-m1, j1+m2-j3) );
kmax = min( j1+j2-j3, min(j1-m1, j2+m2) );
sgn = (-1)^(j1-j2-m3+kmin);

for k = kmin:kmax
   denom = fact(k+1) * fact(j1+j2-j3-k+1) * fact(j1-m1-k+1) * ...
           fact(j2+m2-k+1) * fact(j3-j2+m1+k+1) * fact(j3-j1-m2+k+1);
   result = result + sgn / denom;
   sgn = -sgn;
end
j = j1+j2+j3;
pref = fact(j-2*j1+1) * fact(j-2*j2+1) * fact(j-2*j3+1) / fact(j+2) * ...
       fact(j1+m1+1) * fact(j1-m1+1) * ...
       fact(j2+m2+1) * fact(j2-m2+1) * ...
       fact(j3+m3+1) * fact(j3-m3+1);
result = sqrt( pref ) * result;
end
