
clear all; close all; clc; 

fa = 8e3;
fp = 1050;
fs = 1400;
ws = fs/(fa/2);
wp = fp/(fa/2);
wc = ((wp + ws)/2)*pi;
Ar = 30;
Ap = 0.5;
M = 30;

n = [-M:M];

c_lp = sin(wc*n)./(pi*n);
c_lp(M+1) = wc/pi;

w =  hamming(2*M+1)';
%w = barthannwin(2*M+1)';
%w =  bartlett(2*M+1)';
%w =  hann(2*M+1)';
%w =  rectwin(2*M+1)';
%w =  triang(2*M+1)';
%w =  blackman(2*M+1)';
%w =  blackmanharris(2*M+1)';
%w =  nuttallwin(2*M+1)';
%w =  parzenwin(2*M+1)';

num  = c_lp.*w;

[h,w] = freqz(num,1,20000);
figure(1)
subplot(3,2,1:2)
plot((w/pi)*(fa/2),20*log10(abs(h)))
hold on
grid on
title(['Passa Baixas: com janela de Hamming. Ordem N = ' num2str(M*2)])
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
ylim([-Ar-10 2])
plot([0 fs fs fs],[Ap/2 Ap/2 -Ar -Ar], '--r')
plot([0 fp fp ],-[Ap/2 Ap/2 Ar], '--m')

subplot(3,2,3)
plot((w/pi)*(fa/2),20*log10(abs(h)))
hold on
grid on
title(['Passa Baixas: Atenuação na banda de rejeição'])
axis([1280 1480 -35 -20])
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
plot([0 fs fs fs],[Ap/2 Ap/2 -Ar -Ar], '--r')
plot([0 fp fp ],-[Ap/2 Ap/2 Ar], '--m')

subplot(3,2,4)
plot((w/pi)*(fa/2),20*log10(abs(h)))
hold on
grid on
axis([1049.5 1050.3 -1 0.5])
title(['Passa Baixas: Atenuação na banda de passagem'])
%ylim([-Ar-10 2])
plot([0 fs fs fs],[Ap/2 Ap/2 -Ar -Ar], '--r')
plot([0 fp fp ],-[Ap/2 Ap/2 Ar], '--m')
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');

subplot(3,2,5)
%freqz(num,1);
plot(w*fa/2/pi, angle(h)); grid on;
ylabel('Fase'); xlabel('Frequencia (Hz)');
title('Resposta do filtro - Fase');

subplot(3,2,6)
zplane(num,1);
title('Passa Baixas: Pólos e zeros');
