% DCS for free atoms and for muffin-tin potentials
% see read_eigs.m for more detail regarding parameters

clear
close
tic
load('parameters.mat');
[rad1,dens1,pot1] = denspot(22,90.0);
[rad2,dens2,pot2] = denspot(34,90.0);
rad = {rad1,rad2,rad2};
pot = {pot1,pot2,pot2};
d_ti_se = norm(tau(:,2)-tau(:,1))*alat;
rad1_muff = 0.59*d_ti_se;
rad2_muff = 0.41*d_ti_se;
pot1_muff = muffin(1,2,rad1_muff,rad,pot,tau,at,alat);
pot2_muff = muffin(2,1,rad2_muff,rad,pot,tau,at,alat);

[amp1,~,theta] = scat(90.0,15,rad1,pot1); dcs1 = abs(amp1).^2;
amp2 = scat(90.0,15,rad2,pot2); dcs2 = abs(amp2).^2;
amp1_muff = scat(90.0,15,rad1,pot1_muff);
amp2_muff = scat(90.0,15,rad2,pot2_muff);
dcs1_muff = abs(amp1_muff).^2;
dcs2_muff = abs(amp2_muff).^2;
semilogy(theta,dcs1,'b:',theta,dcs2,'r:','linewidth',2); hold on
semilogy(theta,dcs1_muff,'b-',theta,dcs2_muff,'r-');
legend('free atom Ti','free atom Se','muffin-tin Ti','muffin-tin Se');
xlabel('scattering angle (deg)');
ylabel('DCS (arb)');
set(gca,'FontSize',16);
xlim([0,180]);
grid on

toc