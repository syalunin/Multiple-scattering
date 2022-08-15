function pot_muff = muffin(j,jj,rad_muff,rad,pot,tau,at,alat)
% Generate a muffin-tin potential centred at atom j 
% in the direction of atom jj

rj = tau(:,j)*alat;
rjj = tau(:,jj)*alat;
nj = (rjj-rj)/norm(rjj-rj);

x1 = rj(1)+nj(1)*rad{j};
x2 = rj(2)+nj(2)*rad{j};
x3 = rj(3)+nj(3)*rad{j};

pot_muff = 0*pot{j};
for n1=-7:7
for n2=-7:7
for n3=-7:7
  R = at(:,1)*n1+at(:,2)*n2+at(:,3)*n3;
  for na=1:size(tau,2)
    dx1 = x1-(tau(1,na)+R(1))*alat;
    dx2 = x2-(tau(2,na)+R(2))*alat;
    dx3 = x3-(tau(3,na)+R(3))*alat;
    dist = sqrt(dx1.^2+dx2.^2+dx3.^2);
    pot_muff = pot_muff+spline(rad{na},pot{na},dist);
  end
end
end
end
pot_muff = pot_muff.*(rad{j}<=rad_muff);
end
