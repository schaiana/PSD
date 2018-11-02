close all; clear all; clc;

%1336	1229,135	1313,96	1358,04	1452,845

%% Parametros iniciais

%Atenuacoes de stopband e de passagem (dB)
As = 30;
Ap = 0.5;

%Frequencias de stopband e de passagem (Hz)
fs1 = 1229.135;
fp1 = 1313.96;
fp2 = 1358.04;
fs2 = 1452.845;


%Frequencia de amostragem (Hz)
fa = 3675;

%% Obtencao de parametros do filtro atraves da funcao kaiserord

%ordem minima: n - (fp-fs)/(f_sobra_embaixo + f_sobra_cima)
%verifica se passou alguma coisa e desconta do outro lado (d_f)

d_f = 1; %ajuste de frequencia -> verificado apos alterar a ordem
%f = [fs fp] - d_f;
%mags = [0 1]; %magnitude nas frequencias especificadas em f
f = [fs1 fp1 fp2 fs2] - d_f;
mags = [0 1 0];

devs = [10^(-As/20) (10^(Ap/20)-1)/(10^(Ap/20)+1) 10^(-As/20)];
% devs = [0.1 0.03]; % magnitude do ripple na stopband e na banda de passagem

%beta controla a janela
%Ordem, bandas de freq normalizadas, fator beta e tipo do filtro
[n, Wn, beta, ftype] = kaiserord(f, mags, devs, fa);

%% Projeto do filtro com a funcao fir
n_ajust =  75;%ceil(n - ((1236-1229) / (48+6))); %ajuste da ordem de acordo com a folga analisada no grafico
h = fir1(n_ajust, Wn, ftype, kaiser(n_ajust+1, beta), 'noscale');
[Hw, w] = freqz(h, 1, 2e6);

%% Plots
figure(1)
subplot(3,2,1:2)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([fs1 fs1 fs2 fs2], [-As Ap Ap -As], '--r');
plot([fp1 fp1 fp2 fp2], [-As-20 -Ap -Ap -As-20], '--r');
xlim([1150 1500]); ylim([-60 1]);
%axis([0 3000 -50 2])
%title(['Filtro HP com janela de Kaiser e ordem = ' num2str(n_ajust)]);
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,3)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([fs1 fs1 fs2 fs2], [-As Ap Ap -As], '--r');
plot([fp1 fp1 fp2 fp2], [-As-20 -Ap -Ap -As-20], '--r');
xlim([1305 1365]); ylim([-5 1]);
%axis([1000 1100 -25 -15])
%title('Atenuacao na banda de rejeicao');
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,4)
plot(w*fa/2/pi, 20*log10(abs(Hw)));
grid on; hold on;
plot([fs1 fs1 fs2 fs2], [-As Ap Ap -As], '--r');
plot([fp1 fp1 fp2 fp2], [-As-20 -Ap -Ap -As-20], '--r');
xlim([1200 1475]); ylim([-45 -20]);
%axis([1000 1100 -25 -15])
% axis([1300 1500 -1.5 1.5])
% title('Atenuacao na banda de rejeicao');
ylabel('Magnitude (dB)'); xlabel('Frequencia (Hz)');
hold off;

subplot(3,2,5)
plot(w*fa/2/pi, angle(Hw)); grid on;
ylabel('Fase'); xlabel('Frequencia (Hz)');
title('Resposta do filtro - Fase');

subplot(3,2,6)
zplane(h, 1);
title('Polos e zeros do filtro'); grid on;