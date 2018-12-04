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


%Solucion exacta en el tiempo 0
phi = @(x) cdf('Normal',x,0,1);
Vex = @(S) S.*phi((log(S)+0.14)/sqrt(0.08)) - exp(-0.1)*phi((log(S)+0.06)/sqrt(0.08));

M = @(h) 4/h;
N = @(k) 2/k;

pasos = [1/10 1/20 1/40];
for i = 1:3
   Imp =  mBS_imp(Ix, It, M(pasos(i)), N(pasos(i)), bs);
   gridSpace = linspace(Ix(1), Ix(2), M(pasos(i)) + 1);
   gridPrice = 2.^gridSpace;
   solExact = Vex(gridPrice);
   tamPaso = pasos(i)
   errorMax = max(abs(Imp(:,N(pasos(i))+1)-solExact'))
end


