clear;

load state.dat;
time = state(:,1);
hamiltonian = state(:,2);
potential = state(:,3);
pressure = state(:,4).*16.6;
temperature = state(:,5);

nbins = 100;
m = 39.9;
k = 8.31e-3;
meanT = mean(temperature);

load initial.dat;
pxi = initial(:,4);
pyi = initial(:,5);
pzi = initial(:,6);
momentum = sqrt(pxi.*pxi+pyi.*pyi+pzi.*pzi);

initialT = 2/3/k/size(momentum)(1)/2/m*sum(momentum.*momentum);

load final.dat;
pxf = final(:,4);
pyf = final(:,5);
pzf = final(:,6);
pf = sqrt(pxf.*pxf+pyf.*pyf+pzf.*pzf);

maximumMomentum = max([momentum ; pf]);

momentumRange = 0:0.05:maximumMomentum;
momentumRange2 = power(momentumRange,2);
mkti = 2*m*k*initialT;
mktf = 2*m*k*meanT;
maxwellI = momentumRange2.*exp(-momentumRange2./mkti)*4/power(mkti,1.5)/sqrt(pi).*size(pf)(1)*maximumMomentum/nbins;
maxwellF = momentumRange2.*exp(-momentumRange2./mktf)*4/power(mktf,1.5)/sqrt(pi).*size(pf)(1)*maximumMomentum/nbins;
maximumMomentum

figure(1);
plot(time, hamiltonian, '-', time, potential,'.');
xlabel('Czas [ps]');
ylabel('Energia [kJ/mol]');
legend('Energia calkowita', 'Potencjal');
print('E.png','-S640,500');
close(1);

figure(2);
plot(time, pressure);
xlabel('Czas [ps]');
ylabel('Cisnienie [atm]');
print('P.png','-S640,500');
close(2);

figure(3);
plot(time, temperature);
xlabel('Czas [ps]');
ylabel('Temperatura [K]');
print('T.png','-S640,500');
close(3);

figure(4);
plot(pressure(1),temperature(1),'x0',pressure(end),temperature(end),'o0',pressure, temperature,'-0');
hold on;
xlabel('Cisnienie [atm]');
ylabel('Temperatura [K]');
legend('poczatek','koniec');
hold off;
print('PT.png','-S640,500');
close(4);

final = [0.4, 0.4, 0.4];
initial = [0.8, 0.8, 0.8];

figure(5);
hist(momentum, nbins,'FaceColor',initial,'EdgeColor',initial);
hold on;
hist(pf, nbins,'FaceColor',final,'EdgeColor',final);
plot(momentumRange,maxwellI,'-0',momentumRange,maxwellF,'-0');
legend('ped poczatkowy','ped koncowy', 'rozklad Maxwella');
ylabel('Ilosc czastek');
xlabel('Ped [u nm ps^{-1}]');
print('p.png','-S640,500');
close(5);

disp(['av H: ' num2str(mean(hamiltonian))])
disp(['av P: ' num2str(mean(pressure))])
disp(['av T: ' num2str(meanT)])
