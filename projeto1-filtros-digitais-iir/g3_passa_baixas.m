close all; clear all; clc;

%% Parametros iniciais
fp = 1050; %freq de passagem (Hz)
fs = 1400; %freq de stopband (Hz)
fa = 8e3; %freq de amostragem (Hz)

%Frequencias em radianos
wp = 2*pi*fp;
ws = 2*pi*fs;

G0 = 0; %Ganho na passagem (dB)
Ap = 0.5; %Atenuacao na passagem(dB)
As = 30; %Atenuacao na stopband(dB)

%% Normalizando as frequencias em relacao a freq  de amostragem
% Pre distorcao para compensar a transformada linear

os = ws/fa; Ls = 2*tan(os/2);
op = wp/fa; Lp = 2*tan(op/2);

%% Projeto do passa-baixa Chebychev tipo 1 - Prototipo

%Frequencias de passagem e stopbando do filtro prototipo
W_p = 1;  %wp/wp
W_s = Ls/Lp;

E = sqrt((10^(Ap/10)) - 1); % Define a atenuacao na banda de passagem
n = ceil(acosh(sqrt(((10^(As/10)) - 1) / (E^2))) / (acosh(W_s))); % Ordem 

%Calulo dos polos
k = 1:n;
p = -sinh((1/n)*asinh(1/E))*sin(((2*k - 1)*pi)/(2*n)) + 1j*cosh((1/n)*asinh(1/E))*cos(((2*k - 1)*pi)/(2*n));

if(mod(n,2) == 1)
    G0 = 1;
else
    G0 = sqrt(1/(1 + E^2));
end

den = real(poly(p)); %polos
num = G0*den(end); %zeros

[h, w] = freqs(num, den, 2000000);

%% Transformacao em frequencia do filtro H(p) -> H(s)
% Calculo simbolico
syms s p

%passando numerador e denominador para simbolico
num_simb_p = poly2sym(num, p);
den_simb_p = poly2sym(den, p);

f_s = s / Lp;
%ajusta precisao de f_s para 4 digitos significativos depois da virgula
pretty(vpa(f_s, 4));

Hp(p) = num_simb_p/den_simb_p;
pretty(vpa(Hp(p), 4));

%substitui por f_s todas as ocorrencias de Hp(p)
Hs(s) = subs(Hp(p), f_s); 
pretty(vpa(Hs(s), 4));

%simplifica (enxuga) a equacao
Hsc(s) = collect(Hs(s));
pretty(vpa(Hsc(s), 4));

%extraindo numerador e denominador
[num_sim_s, den_sim_s] = numden(Hsc(s));

%voltando para polinomial e normalizando valores
num_s = sym2poly(num_sim_s);
den_s = sym2poly(den_sim_s);
num_norm = num_s/num_s(1);
den_norm = den_s/num_s(1);

[h_norm, w_norm] = freqs(num_norm, den_norm, 2000000);

%% Transformacao bilinear H(s) -> H(z)

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
Nz_norm = N_z/D_z(1);
Dz_norm = D_z/D_z(1);
num_z_norm = poly2sym(Nz_norm, z);
den_z_norm = poly2sym(Dz_norm, z);
Hzn(z) = num_z_norm/den_z_norm;
pretty(vpa(Hzn(z), 4));

%% Plots - Prototipo

figure(1)
subplot(3,2,6)
zplane(num, den);
title('Pólos e zeros do filtro prototipo H(s)');

subplot(3,2,1:2)
plot(w, 20*log10(abs(h))); xlim([0.8 1.4]); ylim([-40 1]); %modulo de H(p)%Modulo de H(p)
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
hold off
title('Resposta de H(p)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
legend(['Ordem: ' num2str(n)]);

subplot(3,2,3)
plot(w, 20*log10(abs(h))); xlim([0.85 1.05]); ylim([-1 0.2]);
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,4)
plot(w, 20*log10(abs(h))); xlim([1.1 1.7]); ylim([-35 -15]);
hold on; grid on;
plot([0  W_s W_s w(end)], [0 0 -As -As], 'r--');
plot([0 W_p W_p], -[Ap Ap As], 'r--');
title('Atenuação na banda rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,5)
plot(w, angle(h)); xlim([0 2]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% plots - H(s)

figure(2)
subplot(3,2,6)
grid on
zplane(num_norm, den_norm); title ('Pólos e zeros de H(s)');

subplot(3,2,1:2)
plot(w_norm, 20*log10(abs(h_norm))); xlim([0.6 1.3]); ylim([-35 5]);
hold on; grid on;
plot([0 Lp Lp], -[Ap Ap As], 'r--');
plot([0 Ls Ls w_norm(end)],  [0 0 -As -As], 'r--');
title('Resposta de H(s)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off

subplot(3,2,3)
plot(w_norm, 20*log10(abs(h_norm))); ylim([-4 1]); xlim([0.72 0.92]);
hold on; grid on;
plot([0 Lp Lp], -[Ap Ap As], 'r--');
plot([0 Ls Ls w_norm(end)],  [0 0 -As -As], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off

subplot(3,2,4)
plot(w_norm, 20*log10(abs(h_norm))); ylim([-36 -14]); xlim([0.9 1.5]);
hold on; grid on;
plot([0 Lp Lp], -[Ap Ap As], 'r--');
plot([0 Ls Ls w_norm(end)],  [0 0 -As -As], 'r--');
title('Atenuação na banda de rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off

subplot(3,2,5)
plot(w, angle(h)); xlim([0 2]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% Plots H(z)

figure(3)
subplot(3,2,6)
zplane(Nz_norm, Dz_norm); xlim([-2 2]);
title('Pólos e zeros de H(z)');

[hz_norm, wz_norm] = freqz(Nz_norm, Dz_norm, 2000000);

subplot(3,2,1:2)
plot((wz_norm/pi)*(fa/2), 20*log10(abs(hz_norm))); xlim([800 1500]); ylim([-As-10 10]);
hold on; grid on;
plot([0 fp fp], -[Ap Ap As], 'r--');
plot([0 fs fs 1500],[0 0 -As -As], 'r--');
title('Resposta de H(z)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,3)
plot((wz_norm/pi)*(fa/2), 20*log10(abs(hz_norm))); xlim([900 1100]); ylim([-2 0.1]);
hold on; grid on;
plot([0 fp fp], -[Ap Ap As], 'r--');
plot([0 fs fs 1500],[0 0 -As -As], 'r--');
title('Atenuação na banda de passagem'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,4)
plot((wz_norm/pi)*(fa/2), 20*log10(abs(hz_norm))); xlim([1360 1440]); ylim([-40 -28]);
hold on; grid on;
plot([0 fp fp], -[Ap Ap As], 'r--');
plot([0 fs fs 1500],[0 0 -As -As], 'r--');
title('Atenuação na banda de rejeição'); ylabel('Atenuação (dB)'); xlabel('banda de passagem');
hold off;

subplot(3,2,5)
plot(wz_norm, angle(hz_norm)); xlim([0 1])
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('banda de passagem');

%% Atraso de grupo do Filtro

figure(4)
grpdelay(Nz_norm, Dz_norm)
title('Atraso de grupo de H(z)');