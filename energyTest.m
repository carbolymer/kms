clear;
energy1 = load('-ascii','t1e-3/state.dat')(:,2);
energy1 = energy1./energy1(1);
energy2 = load('-ascii','t2e-3/state.dat')(:,2);
energy2 = energy2./energy2(1);
energy3 = load('-ascii','t5e-4/state.dat')(:,2);
energy3 = energy3./energy3(1);

k = 1:1:size(energy1)(1);
k = k.*5;

figure(1);
hold on;
plot(k,energy2,k,energy1,k,energy3);

hold off;
axis([0 k(end) 1 1.6]);
ylabel('Wzrost energi H_{k}/H_{k=0}');
xlabel('Numer kroku k');
legend('\tau = 2\cdot10^{-3} ps', '\tau = 1\cdot10^{-3} ps', '\tau = 5\cdot10^{-4} ps');
print('stability.png','-S640,500');
close(1);