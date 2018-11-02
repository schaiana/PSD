close all; clear all; clc;

%% Parametros iniciais
fp = 1400; %freq de passagem (Hz)
fs = 1050; %freq de stopband (Hz)
fa = 8e3; %freq de amostragem (Hz)

Ap = 0.5; %Atenuacao na passagem(dB)
As = 20; %Atenuacao na stopband(dB)

%% Normalizando as frequencias em relacao a freq  de amostragem
% Pre distorcao para compensar a transformada linear

os = 2*pi*fs/fa; Ls = 2*tan(os/2);
op = 2*pi*fp/fa; Lp = 2*tan(op/2);

%% Projeto do passa-alta Butterworth - Prototipo

%Frequencias de passagem e stopbando do filtro prototipo
W_p = 1;  %wp/wp
W_s = Lp/Ls;

E = sqrt(10^(Ap/10)-1); % Define a atenuacao na banda de passagem

n = ceil((log((10^(0.1*As) - 1) / (10^(0.1*Ap) - 1))) / (2*log(W_s)));% Ordem 

%Calulo dos polos
k = 1:n;
p = E^(-1/n) * exp(1j*(((2*k)+n-1)/(2*n))*pi);

%Obtendo H(p)
a = real(poly(p));
b = a(end);

syms p
Np = poly2sym(b,p);
Dp = poly2sym(a,p);
Hp(p) = Np/Dp;

[h,w] = freqs(b, a, 2e6);

%% Transformacao em frequencia do filtro H(p) -> H(s) Passa alta
% Calculo simbolico
syms s p

f_s = Lp/s;
W_s = Lp/Ls;

Hs(s) = subs(Hp(p), f_s);
Hsc(s) = collect(Hs(s));
[Ns, Ds] = numden(Hsc(s));

num_s = sym2poly(Ns);
den_s = sym2poly(Ds);

%normalizacao do numerador e denominador
num_norm = num_s/den_s(1);
den_norm = den_s/den_s(1);

[h_s,w_s] = freqs(num_norm, den_norm, 2e6);

%% Transformacao bilinear H(s) -> H(z)

T = 1/fa; %periodo de amostragem

syms z
num_norm_s = poly2sym(num_norm, s);
den_norm_s = poly2sym(den_norm, s);

f_z = 2*(z-1)/(z+1);
pretty(vpa(f_z, 4));
Hs_norm(s) = num_norm_s/den_norm_s;
pretty(vpa(Hs_norm(s), 4));

Hz(z) = subs(Hs_norm(s), f_z);
pretty(vpa(Hz(z), 4));
Hzc(z) = collect(Hz(z));
pretty(vpa(Hzc(z), 4));

[Nz, Dz] = numden(Hzc(z));
N_z = sym2poly(Nz);
D_z = sym2poly(Dz);

%normalizacao do numerador e denominador
Nz_norm = N_z/D_z(1);
Dz_norm = D_z/D_z(1);
num_z_norm = poly2sym(Nz_norm, z);
den_z_norm = poly2sym(Dz_norm, z);
Hzn(z) = num_z_norm/den_z_norm;
pretty(vpa(Hzn(z), 4));

%% Plots - Prototipo H(p)

figure(1)
subplot(3,2,6)
zplane(b, a);
title('Pólos e zeros do filtro prototipo H(p)');

subplot(3,2,1:2)
plot(w, 20*log10(abs(h))); xlim([0.6 1.8]); ylim([-30 1]); %modulo de H(p)
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
hold off
title('Resposta de H(p)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
legend(['Ordem: ' num2str(n)]);

subplot(3,2,3)
plot(w, 20*log10(abs(h))); xlim([0.99 1.01]); ylim([-1 0.2]);
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,4)
plot(w, 20*log10(abs(h))); xlim([1.38 1.43]); ylim([-25 -15]);
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
title('Atenuação na banda rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,5)
plot(w, angle(h)); xlim([0 2]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% Plots - H(s)

figure(2)
subplot(3,2,6)
zplane(num_norm, den_norm);
title('Pólos e zeros do filtro prototipo H(s)');

subplot(3,2,1:2)
plot(w_s, 20*log10(abs(h_s))); xlim([0.6 1.6]); ylim([-30 1]); %modulo de H(s)
hold on; grid on;
plot([0  Ls Ls w_s(end)], [-As -As Ap Ap], 'r--');
plot([Lp Lp w_s(end)], -[30 Ap Ap], 'r--');
hold off
title('Resposta de H(s)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
legend(['Ordem: ' num2str(n)]);

subplot(3,2,3)
plot(w_s, 20*log10(abs(h_s))); xlim([1.2 1.26]); ylim([-0.9 -0.3]); %modulo de H(s)
hold on; grid on;
plot([0  Ls Ls w_s(end)], [-As -As Ap Ap], 'r--');
plot([Lp Lp w_s(end)], -[30 Ap Ap], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,4)
plot(w_s, 20*log10(abs(h_s))); xlim([0.8 0.9]); ylim([-24 -18]); %modulo de H(s)
hold on; grid on;
plot([0  Ls Ls w_s(end)], [-As -As Ap Ap], 'r--');
plot([Lp Lp w_s(end)], -[30 Ap Ap], 'r--');
title('Atenuação na banda rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,5)
plot(w_s, angle(h_s)); xlim([0 4]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% Plots - H(z)

figure(3)
subplot(3,2,6)
zplane(Nz_norm, Dz_norm); xlim([-2 2]);
title('Pólos e zeros de H(z)');

[hz, wz] = freqz(Nz_norm, Dz_norm, 2e6);

subplot(3,2,1:2)
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); xlim([900 1600]); ylim([-25 2]);
hold on; grid on;
plot([0  fs fs 1600], [-As -As Ap Ap], 'r--');
plot([fp fp 1600], -[30 Ap Ap], 'r--');
title('Resposta de H(z)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,3)
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); xlim([1350 1600]); ylim([-4 2]);
hold on; grid on;
plot([0  fs fs 1600], [-As -As Ap Ap], 'r--');
plot([fp fp 1600], -[30 Ap Ap], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,4)
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); xlim([1000 1100]); ylim([-24 -16]);
hold on; grid on;
plot([0  fs fs 1600], [-As -As Ap Ap], 'r--');
plot([fp fp 1600], -[30 Ap Ap], 'r--');
title('Atenuação na banda de rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,5)
plot(wz, angle(hz)); xlim([0 2])
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% Atraso de grupo do Filtro
figure(4)
grpdelay(Nz_norm, Dz_norm)
title('Atraso de grupo de H(z)');