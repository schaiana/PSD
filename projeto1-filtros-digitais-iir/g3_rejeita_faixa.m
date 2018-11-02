close all; clear all; clc;

%% Atenuacoes e Frequencias de delimitacao
AtenPass = 0.5; %dB
AtenStop = 30; %dB

freqPass1 = 767;
freqStop1 = 835; %Hz
freqStop2 = 869;
freqPass2 = 937;

fa = 8e3; %Hz

%Normalizando as freq em relação a fa
os1 = freqStop1*2*pi/fa; Ls1 = 2*tan(os1/2); %pre distorcao para compensar a transf bilinear
os2 = freqStop2*2*pi/fa; Ls2 = 2*tan(os2/2);
op1 = freqPass1*2*pi/fa; Lp1 = 2*tan(op1/2);
op2 = freqPass2*2*pi/fa; Lp2 = 2*tan(op2/2);

%f0 = 852; valor teorico - freq de passagem
L0 = sqrt(Lp1 * Lp2);
w0 = sqrt(freqPass1*freqPass2) * 2*pi;

BL = (Lp2 - Lp1);
%% Calculo do Filtro passa-baixa
Ws1 = abs((Ls1*BL)/((L0^2 - Ls1^2)));
Ws2 = abs((Ls2*BL)/((L0^2 - Ls2^2)));
Wstop1 = min(Ws1, Ws2); %H(p)
Wpass1 = 1;

[n,Wp] = ellipord(Wpass1,Wstop1,AtenPass,AtenStop,'s');
[b,a] = ellip(n,AtenPass,AtenStop,Wp,'s');
[h,w] = freqs(b,a,fa);

%% Plot protótipo H(p)
figure(1)
subplot(3,2,1:2)
plot(w, 20*log10(abs(h))); xlim([-0.1 2]); ylim([-20 1]); legend(['Filtro Eliptico. Ordem:' num2str(n)]);
hold on
plot([Wstop1 Wstop1+5], -[AtenStop AtenStop], 'm--')
plot([0 Wstop1 Wstop1],  [0 0 -AtenStop], 'm--');
plot([0 Wpass1 Wpass1], -[AtenPass AtenPass AtenStop], 'm--');
legend(['Filtro Eliptico: ' num2str(n)]);
title('Resposta de H(p)- Magnitude'); ylabel('Atenuação (dB)'); xlabel('Banda de passagem[Hz]');
ylim([-35 0]);
xlim([0 5]);

subplot(3,2,3)
plot(w, 20*log10(abs(h))); xlim([-0.1 2]); ylim([-20 1]); 
hold on
plot([0 Wstop1 Wstop1],  [0 0 -AtenStop], 'm--');
plot([0 Wpass1 Wpass1], -[AtenPass AtenPass AtenStop], 'm--');
title('Atenuação na banda de passagem'); xlabel('Banda Passante [Hz]'); ylabel('Atenuação [dB]');
ylim([-4 0]);
xlim([0.95 1.05]);

subplot(3,2,4)
plot(w, 20*log10(abs(h))); xlim([-0.1 2]); ylim([-20 1]); 
hold on
plot([0 Wstop1 Wstop1],  [0 0 -AtenStop], 'm--');
plot([Wstop1 Wstop1+5], -[AtenStop AtenStop], 'm--')
plot([0 Wpass1 Wpass1], -[AtenPass AtenPass AtenStop], 'm--');
title('Atenuação na banda rejeição'); ylabel('Atenuação [dB]'); xlabel('Banda de passagem [Hz]');
ylim([-35 -20]);
xlim([Wstop1-0.05 Wstop1+0.05]);

subplot(3,2,5)
plot(w, angle(h)); xlim([-0.1 2]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('Banda de passagem [Hz]');

subplot(3,2,6)
zplane(b,a);
title('Polos e zeros do protótipo H(p)');

%% C?lculo H(s)
syms s p
Np = poly2sym(b, p);
Dp = poly2sym(a, p);

f_s = (BL*s)/(s^2 + L0^2);


%1pretty(vpa(f_s, 4));
Hp(p) = Np/Dp;
pretty(vpa(Hp(p), 4));
Hs(s) = subs(Hp(p), f_s);
%3pretty(vpa(Hs(s), 4));

Hsc(s) = collect(Hs(s));
%4pretty(vpa(Hsc(s), 4));

[Ns, Ds] = numden(Hsc(s));
bs = sym2poly(Ns);
as = sym2poly(Ds);

as_norm = as/as(1);
bs_norm = bs/as(1);

figure(2)
subplot(3,2,1:2); hold on; grid on;
[h, w] = freqs(bs_norm, as_norm, 20000); 
plot(w, 20*log10(abs(h)));
plot([Lp1-1 Lp1 Lp1], -[AtenPass AtenPass AtenStop], '--r')
plot([Lp1-1 Ls1 Ls1 Ls2 Ls2 Lp2+1], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([Lp2+1 Lp2 Lp2], -[AtenPass AtenPass AtenStop], '--r')
title('Filtro analógico rejeita faixa - Protótipo H(s)'); xlabel('Banda Passante [Hz]'); ylabel('Magnitude [dB]');
xlim([Lp1-0.01 Lp2+0.01]);
ylim(-[AtenStop+6 AtenPass-5]);
hold off;

subplot(3,2,3);
plot(w, 20*log10(abs(h)));
hold on; grid on;
plot([Lp1-1 Lp1 Lp1], -[AtenPass AtenPass AtenStop], '--r')
plot([Lp1-1 Ls1 Ls1 Ls2 Ls2 Lp2+1], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([Lp2+1 Lp2 Lp2], -[AtenPass AtenPass AtenStop], '--r')
title('Atenuação na banda de passagem'); xlabel('Banda Passante [Hz]'); ylabel('Atenuação [dB]');xlim([Lp1-0.05 Lp2+0.05]);
ylim(-[AtenPass+1 AtenPass-1]);
hold off;

subplot(3,2,4);
plot(w, 20*log10(abs(h)));
hold on; grid on;
plot([Lp1-1 Lp1 Lp1], -[AtenPass AtenPass AtenStop], '--r')
plot([Lp1-1 Ls1 Ls1 Ls2 Ls2 Lp2+1], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([Lp2+1 Lp2 Lp2], -[AtenPass AtenPass AtenStop], '--r')
title('Atenuação na banda rejeição'); ylabel('Atenuação [dB]'); xlabel('Banda de passagem [Hz]');
xlim([Lp1-0.05 Lp2+0.05]);
ylim(-[AtenStop+1 AtenStop-1]);
hold off;


subplot(3,2,5)
plot(w, angle(h)); %xlim([-0.1 2]);
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('Banda de passagem [Hz]');

subplot(3,2,6)
zplane(bs_norm,as_norm);
title('Polos e zeros do protótipo H(p)');

%% Transforma??o bilinear
T = 1;

syms z
Ns_norm = poly2sym(bs_norm,s);
Ds_norm = poly2sym(as_norm,s);
f_z = 2*(z - 1)/(T*(z+1));
Hs_norm(s) = Ns_norm/Ds_norm;
pretty(vpa(Hs_norm(s), 4));

Hz(z) = subs(Hs_norm(s), f_z);
%6pretty(vpa(Hz(z), 4));
Hzc(z) = collect(Hz(z));
%pretty(vpa(Hzc(z), 4));

%% H(z)
[Nz, Dz] = numden(Hzc(z));
bz = sym2poly(Nz);
az = sym2poly(Dz);
azn = az/az(1);
bzn = bz/az(1);
Nzn = poly2sym(bzn,z);
Dzn = poly2sym(azn,z);
Hzn(z) = Nzn/Dzn;
%pretty(vpa(Hzn(z), 4));

figure(3)
subplot(3,2,1:2);
[hz, wz] = freqz(bzn, azn,15000);
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); 
hold on; grid on;
plot([freqPass1-100 freqPass1 freqPass1], -[AtenPass AtenPass AtenStop], '--r')
plot([freqPass1-100 freqStop1 freqStop1 freqStop2 freqStop2 freqPass2+100], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([freqPass2+100 freqPass2 freqPass2], -[AtenPass AtenPass AtenStop], '--r')
xlim([freqPass1-50 freqPass2+50]);
ylim(-[AtenStop+5 AtenPass-5])
title('Resposta de H(z)- Magnitude')
xlabel('Banda Passante [Hz]'); ylabel('Atenuação [dB]');
hold off;

subplot(3,2,3);
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); 
hold on; grid on;
plot([freqPass1-100 freqPass1 freqPass1], -[AtenPass AtenPass AtenStop], '--r')
plot([freqPass1-100 freqStop1 freqStop1 freqStop2 freqStop2 freqPass2+100], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([freqPass2+100 freqPass2 freqPass2], -[AtenPass AtenPass AtenStop], '--r')
xlim([freqPass1-10 freqPass2+10]); ylim(-[AtenStop+1 AtenStop-1]);
title('Atenuação na banda de passagem'); xlabel('Banda Passante [Hz]'); ylabel('Atenuação [dB]');

hold off;

subplot(3,2,4);
[hz, wz] = freqz(bzn, azn,15000);
plot((wz/pi)*(fa/2), 20*log10(abs(hz))); 
hold on; grid on;
plot([freqPass1-100 freqPass1 freqPass1], -[AtenPass AtenPass AtenStop], '--r')
plot([freqPass1-100 freqStop1 freqStop1 freqStop2 freqStop2 freqPass2+100], [0 0 -AtenStop -AtenStop 0 0], '--m')
plot([freqPass2+100 freqPass2 freqPass2], -[AtenPass AtenPass AtenStop], '--r')
xlim([freqPass1-10 freqPass2+10]); ylim(-[AtenPass+1 AtenPass-1]);
title('Atenuação na banda rejeição'); ylabel('Atenuação [dB]'); xlabel('Banda de passagem [Hz]');
xlim([freqPass1-50 freqPass2+50]); ylim(-[AtenPass+1 AtenPass-1]);
hold off;

subplot(3,2,5)
plot(wz, angle(hz)); xlim([0 2])
grid on;
title ('Resposta do filtro - Fase'); ylabel('Variação da Fase'); xlabel('Banda de passagem[dB]');


subplot(3,2,6)
zplane(bzn,azn);
grid on;
title('Polos e zeros de H(z)');

figure(4)
grpdelay(bzn, azn)
title('Atraso de grupo de H(z)');
