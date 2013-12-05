clear;
L = 2.3;
coefficient = 3/2*125*8.31e-3/(4/3*pi*power(L,3)).*16.6;
T = [137.766 519.67 1071.21 1582.79];
P = [20.019 324.424 542.231 687.848];
T_teo = 100:200:1800;
P_teo = T_teo.*coefficient;
figure(1);
plot(T,P,'*;Dane z symulacji;',T_teo,P_teo,'-;Teoria;')
ylabel('Srednie cisnienie [atm]');
xlabel('Srednia temperatura [K]');
print('pvnkt.png','-S640,500');
close(1);


figure(2);
a = [36 36.5 37 37.25 37.5 37.75 38 39];
H = [-635.9 -662.3 -674.9 -676.3 -675.8 -673.3 -669 -638.8];
plot(a,H);
ylabel('H [kJ / mol]');
xlabel('a [nm]');
print('net_const.png','-S640,500');
close(2);