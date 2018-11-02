close all; clear all; clc;

%% Parametros iniciais

%Atenuacoes de stopband e de passagem (dB)
As = 20;
Ap = 0.5;

%Frequencias de stopband e de passagem (Hz)
fs = 1050;
fp = 1400;

%Frequencia de amostragem (Hz)
fa = 8e3;

%% Obtencao de parametros do filtro atraves da funcao kaiserord

%ordem minima: n - (fp-fs)/(f_sobra_embaixo + f_sobra_cima)
%verifica se passou alguma coisa e desconta do outro lado (d_f)

d_f = 25; %ajuste de frequencia -> verificado apos alterar a ordem
f = [fs fp] - d_f;
mags = [0 1]; %magnitude nas frequencias especificadas em f
devs = [10^(-As/20) (10^(Ap/20)-1)/(10^(Ap/20)+1)];
% devs = [0.1 0.03]; % magnitude do ripple na stopband e na banda de passagem

%beta controla a janela
%Ordem, bandas de freq normalizadas, fator beta e tipo do filtro
[n, Wn, beta, ftype] = kaiserord(f, mags, devs, fa);

%% Projeto do filtro com a funcao fir
n_ajust = ceil(n - ((fp -fs) / (48+6))); %ajuste da ordem de acordo com a folga analisada no grafico
h = fir1(n_ajust, Wn, ftype, kaiser(n_ajust+1, beta), 'noscale');
[Hw, w] = freqz(h, 1, 2e3);

%% Plots
figure(1)
subplot(3,2,1:2)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([0 fs fs 3000], [-As -As Ap/2 Ap/2], '--r');
plot([fp fp 3000], -[50 Ap/2 Ap/2], '--r');
axis([0 3000 -50 2])
title(['Filtro HP com janela de Kaiser e ordem = ' num2str(n_ajust)]);
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,3)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([0 fs fs 3000], [-As -As Ap/2 Ap/2], '--r');
plot([fp fp 3000], -[50 Ap/2 Ap/2], '--r');
axis([1000 1100 -25 -15])
title('Atenuacao na banda de rejeicao');
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,4)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([0 fs fs 3000], [-As -As Ap/2 Ap/2], '--r');
plot([fp fp 3000], -[50 Ap/2 Ap/2], '--r');
axis([1300 1500 -1.5 1.5])
title('Atenuacao na banda de rejeicao');
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,5)
plot(w*fa/2/pi, angle(Hw)); grid on;
ylabel('Fase'); xlabel('Frequencia (Hz)');
title('Resposta do filtro - Fase');

subplot(3,2,6)
zplane(h, 1);
title('Polos e zeros do filtro'); grid on;