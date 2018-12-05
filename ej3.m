%EJERCICIO 2
close all; clear all; clc;

%Grafica de V(S,t)
Ix = [-2 2];
It = [0 2];
M = 80;
N = 20;

bs.r = 0.05;
bs.sigma = 0.2;
bs.fc = @(x) max([2^x-1 0]);
bs.bcL = @(t) 0;
bs.bcR = @(t) 4*log(2);

W = mBS_imp(Ix, It, M, N, bs);
gridSpace = linspace(Ix(1), Ix(2), M + 1);
gridTime  = linspace(It(1), It(2), N + 1);
gridPrice = 2.^gridSpace;

mesh(gridPrice, gridTime, W', 'LineWidth', 1.5);

%%
%Solucion exacta en el tiempo 0
%phi = @(x) cdf('Normal',x,0,1);
phi = @(x) normcdf(x);
Vex = @(S) S.*phi((log(S)+0.14)/sqrt(0.08)) - exp(-0.1)*phi((log(S)+0.06)/sqrt(0.08));

M = @(h) 4/h;
N = @(k) 2/k;

pasos = [1/10 1/20 1/40];
for i = 1:3
   Imp =  suave(Ix, It, M(pasos(i)), N(pasos(i)), bs);
   gridSpace = linspace(Ix(1), Ix(2), M(pasos(i)) + 1);
   gridPrice = 2.^gridSpace;
   solExact = Vex(gridPrice);
   tamPaso = pasos(i)
   errorMax = max(abs(Imp(:,N(pasos(i))+1)-solExact'))
end
%pause();
%Inciso b)
%Utilizando Crank-Nicolson graficamos las soluciones analiticas y numericas

Crank =  suave(Ix, It, M(1/100), N(1/10), bs);
tFin = N(1/10)+1;
gridSpace = linspace(Ix(1), Ix(2), M(1/100) + 1);
gridPrice = 2.^gridSpace;
solExact = Vex(gridPrice);

%Grafica de s contra V numerica y analitica
close all;
%plot(gridPrice,solExact)
hold on
%plot(gridPrice,Crank(:,N(1/10)+1))
%plot(gridPrice,Crank(:,1))
%pause();

%Grafica de Vs contra s numerica y analitica
dn = @(z) 1/sqrt(2*pi)*exp((-z.^2)/2);
Vs = @(s) phi((log(s)+0.14)/sqrt(0.08))+dn((log(s)+0.14)/sqrt(0.08))/sqrt(0.08)-exp(-0.1)*dn((log(s)+0.06)/sqrt(0.08))/(s*sqrt(0.08));
%VS numérica
%Wx_anal = @(x) Vs(2.^x).*2.^x*log(2);
%Calculamos Wx num�rica
Wx_num(1) = 0;
for i = 2:M(1/100)
   Wx_num(i) = (Crank(i+1,tFin)-Crank(i-1,tFin))*5;
end

solExact = Vex(linspace(0,4,M(1/100)+1));
for i = 2:M(1/100)
   Vs_anal(i) = (solExact(i+1)-solExact(i-1))*5;
end

Vs_anal(1)=0;
Vs_anal(M(1/100)+1) = 1;

Wx_num(M(1/100)+1) = 2.8;
Vs_num = @(Wx,x) Wx./(2.^(x)*log(2));
Vs_res = Vs_num(Wx_num,gridSpace);
close all;
title('Vs vs S')
plot(gridPrice,Vs_res)
hold on
plot(linspace(0,4,M(1/100)+1),Vs_anal)

%Vxx analitica y numerica
Wxx_num(1) = 0;
for i = 2:M(1/100)
   Wxx_num(i) = (Wx_num(i+1)-Wx_num(i-1))*5;
end

for i = 2:M(1/100)
   Vss_anal(i) = (Vs_anal(i+1)-Vs_anal(i-1))*5;
end

Wxx_num(M(1/100)+1) = 2.8;
Vss_anal(1)=0;
Vss_anal(M(1/100)+1) = 0;
Vss_num = @(Wx,Wxx,x) Wxx*(1/(2^x*log(2)))^2 + Wx*(-1/((2^x)^2*log(2)));
for j=1:M(1/100)+1
   Vss_res(j) = Vss_num(Wx_num(j),Wxx_num(j),gridSpace(j));
end
%close all;
plot(gridPrice,Vss_res)
hold on
plot(linspace(0,4,M(1/100)+1),Vss_anal)