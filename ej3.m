%EJERCICIO 3
close all; clear all; clc;

%Parámetros
Ix = [-2 2];
It = [0 2];
M = @(h) 4/h;
N = @(k) 2/k;

bs.r = 0.05;
bs.sigma = 0.2;
bs.fc = @(x) max([2^x-1 0]);
bs.bcL = @(t) 0;
bs.bcR = @(t) 4*log(2);
phi = @(x) normcdf(x);
Vex = @(S) S.*phi((log(S)+0.14)/sqrt(0.08)) - exp(-0.1)*phi((log(S)+0.06)/sqrt(0.08));
%Graficas
%Utilizando Suavizamiento graficamos las soluciones analiticas y numericas

solSuave =  suave(Ix, It, M(1/100), N(1/10), bs);
tFin = N(1/10)+1;
gridSpace = linspace(Ix(1), Ix(2), M(1/100) + 1);
gridPrice = 2.^gridSpace;
solExact = Vex(gridPrice);

%Grafica de s contra V numerica y analitica
close all;
plot(gridPrice,solExact)
hold on
plot(gridPrice,solSuave(:,N(1/10)+1))
plot(gridPrice,solSuave(:,1))
pause();

%Grafica de Vs contra s numerica y analitica
dn = @(z) 1/sqrt(2*pi)*exp((-z.^2)/2);
Vs = @(s) phi((log(s)+0.14)/sqrt(0.08))+dn((log(s)+0.14)/sqrt(0.08))/sqrt(0.08)-exp(-0.1)*dn((log(s)+0.06)/sqrt(0.08))/(s*sqrt(0.08));

Wx_num(1) = (-3*solSuave(1,tFin)+4*solSuave(2,tFin)-solSuave(3,tFin))*50
Wx_num(M(1/100)+1) = (-solSuave(M(1/100)-1,tFin)+4*solSuave(M(1/100),tFin)-3*solSuave(M(1/100)+1,tFin))*-50;

for i = 2:M(1/100)
   Wx_num(i) = (solSuave(i+1,tFin)-solSuave(i-1,tFin))*50;
end

solExact = Vex(linspace(0,4,M(1/100)+1));
for i = 2:M(1/100)
   Vs_anal(i) = (solExact(i+1)-solExact(i-1))*50;
end

%Vs_anal(1)=0;
%Vs_anal(M(1/100)+1) = 1;
Vs_anal(1) = (-3*solExact(1)+4*solExact(2)-solExact(2))*50;
Vs_anal(M(1/100)+1) = (-solExact(M(1/100)-1)+4*solExact(M(1/100))-3*solExact(M(1/100)+1))*-50;


Vs_num = @(Wx,x) Wx./(2.^(x)*log(2));
Vs_res = Vs_num(Wx_num,gridSpace);
close all;
title('Vs vs S')
plot(gridPrice,Vs_res)
hold on
plot(linspace(0,4,M(1/100)+1),Vs_anal)
pause();

%Vxx analitica y numerica
Wxx_num(1) = 0;
for i = 2:M(1/100)
   Wxx_num(i) = (Wx_num(i+1)-Wx_num(i-1))*50;
end

for i = 2:M(1/100)
   Vss_anal(i) = (Vs_anal(i+1)-Vs_anal(i-1))*50;
end

Wxx_num(M(1/100)+1) = 2.8;
Vss_anal(1)=0;
Vss_anal(M(1/100)+1) = 0;
Vss_num = @(Wx,Wxx,x) Wxx*(1/(2^x*log(2)))^2 + Wx*(-1/((2^x)^2*log(2)));
for j=1:M(1/100)+1
   Vss_res(j) = Vss_num(Wx_num(j),Wxx_num(j),gridSpace(j));
end
close all;
plot(gridPrice,Vss_res)
hold on
plot(linspace(0,4,M(1/100)+1),Vss_anal)