% Compute spherical harmonic
% 4-Aug-2022

function result = spherical_harmonic(l,m,theta,phi)
  if abs(m) <= l
    N_lm = legendre(l, cos(theta), 'norm');
    N_lm = reshape(N_lm(abs(m)+1,:), size(theta));
  else
    N_lm = zeros(size(theta));
  end
  result = exp(1i*m*phi) .* N_lm / sqrt(2*pi);
  if m > 0
    result = (-1)^m * result;
  end
end
