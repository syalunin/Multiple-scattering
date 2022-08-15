% plot electron density distribution

clear
close
tic
[rad1,dens1,pot1] = denspot(22,90.0);
[rad2,dens2,pot2] = denspot(34,90.0);

subplot(1,2,1)
semilogx(rad1,rad1.^2.*dens1,'b-'); hold on
semilogy(rad2,rad2.^2.*dens2,'r-');
set(gca,'FontSize',16);
xlabel('rad (au)');
ylabel('density');
legend('Ti','Se');
pbaspect([2 3 1]);
xlim([1e-3,10]);
grid on

subplot(1,2,2)
semilogy(rad1,dens1,'b-',rad2,dens2,'r-');
set(gca,'FontSize',16);
xlabel('rad (au)');
ylabel('density');
legend('Ti','Se');
pbaspect([2 3 1]);
ylim([1e-3,1e3]);
xlim([0,3]);
grid on