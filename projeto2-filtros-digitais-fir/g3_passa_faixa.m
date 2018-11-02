%% Aula 12/05 - BANDPASS with Parks-McClellan (TERMINAR!!!!) 
%BP1 - (f1 = 693 Hz; f2 = 755 Hz, f3 = 785 Hz; f4 = 847 Hz, Elíptico )
close all; clear all; clc;

Ap = 0.5;
As = 30;
fs1 = 693;
fp1 = 755;
fp2 = 785;
fs2 = 847;
rpp = Ap;
rsp = As;
fs = 8000;
f = [fs1 fp1 fp2 fs2];
a = [0 0.99 0];
facDens = 50; %fator de densidade
%magnitude do ripple na stopband e na banda de passagem
dev = [(10^(-As/20))*2.45 (10^(Ap/20)-1)/(10^(Ap/20)+1)*0.1 (10^(-As/20))*2.45];

[n, fo, ao, w] = firpmord(f, a, dev, fs);
b = firpm(n, fo, ao, w);

[hz, wz] = freqz(b, 1, 2048, fs);
figure(1);
subplot(3,2,1:2)
plot(wz, 20*log10(abs(hz))); hold on;
plot([0 fs1 fs1 fs2 fs2 2000], [-As -As 0 0  -As -As], 'm--');
plot([fp1 fp1 fp2 fp2], -[As Ap Ap As], 'r--');
xlim([fs1-100 fs2+100]);
ylim([-33 0.5]);
title(['Parks-McClellan Band-Pass Filter order = ' num2str(length(b)-1)]);

subplot(3,2,3)
plot(wz, 20*log10(abs(hz))); hold on;
plot([0 fs1 fs1 fs2 fs2 2000], [-As -As 0 0  -As -As], 'm--');
plot([fp1 fp1 fp2 fp2], -[As Ap Ap As], 'r--');
xlim([844 850]);
ylim([-32 -27]);
title('Frequência na banda de rejeição');

subplot(3,2,4)
plot(wz, 20*log10(abs(hz))); hold on;
plot([0 fs1 fs1 fs2 fs2 2000], [-As -As 0 0  -As -As], 'm--');
plot([fp1 fp1 fp2 fp2], -[As Ap Ap As], 'r--');
xlim([745 794]);
ylim([-1 0.5]);
title('Frequência na banda de passagem');

subplot(3,2,5)
plot(wz,angle(hz));
title('Resposta do Filtro - fase');
xlim([0 2000]);

subplot(3,2,6)
zplane(b,1);
title('Polos e zeros do Filtro Digital H(z)');
