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
%phi = @(x) cdf('Normal',x,0,1);
phi = @(x) normcdf(x);
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
pause();
%Inciso b)
%Utilizando Crank-Nicolson graficamos las soluciones analiticas y numericas

Crank =  mBS_CN(Ix, It, M(1/10), N(1/10), bs);
tFin = N(1/10)+1;
gridSpace = linspace(Ix(1), Ix(2), M(1/10) + 1);
gridPrice = 2.^gridSpace;
solExact = Vex(gridPrice);

%Grafica de s contra V numerica y analitica
close all;
plot(gridSpace,solExact)
hold on
plot(gridSpace,Crank(:,N(1/10)+1))
plot(gridSpace,Crank(:,1))

%Grafica de Wx contra x numerica y analitica
dn = @(z) 1/sqrt(2*pi)*exp((-z.^2)/2);
Vs = @(s) phi((log(s)+0.14)/sqrt(0.08))+dn((log(s)+0.14)/sqrt(0.08))/sqrt(0.08)-exp(-0.1)*dn((log(s)+0.06)/sqrt(0.08))/(s*sqrt(0.08));
%Wx analitica
Wx_anal = @(x) Vs(2.^x).*2.^x*log(2);
%Calculamos Wx num�rica
Wx_num(1) = 0;
for i = 2:M(1/10)
    Wx_num(i) = (Crank(i+1,tFin)-Crank(i-1,tFin))*5;
end
Wx_num(M(1/10)+1) = 2.8;
%pause();
%close all;
plot(gridSpace,Wx_anal(gridSpace));
hold on
%plot(gridSpace,Wx_num);

%Grafica Wxx contra x numerica y analitica
dn_p = @(z) dn(z).*(-z);
Vss = @(s) dn_p((log(s)+0.14)/sqrt(0.08))./(0.08*s)+dn((log(s)+0.14)/sqrt(0.08))./(sqrt(0.08)*s)-dn_p((log(s)+0.06)/sqrt(0.08))*exp(-0.1)./(s.^2*0.08)+dn_p((log(s)+0.06)/sqrt(0.08))*exp(-0.1)./(s.^2*sqrt(0.08)); 



%dn( (log(s)+0.14)/sqrt(0.08))./(s*sqrt(0.08)) + dn_p((log(s)+0.14)/sqrt(0.08))./(s*0.08)-exp(-0.1)*dn_p((log(s)+0.06)/sqrt(0.08))./(s*0.08);
Wxx_anal = @(x) Vss(2.^x).*(2.^x*log(2)).^2+2.^x*log(2)^2.*Vs(2.^x);
close all;
plot(gridSpace, Wxx_anal(gridSpace));
%pause();
bander = 0
%Prueba Vs y Vss
close all;
p = linspace(1,2,200);
%plot(p,Vs(p))
hold on
%plot(p,Vss(p),'o')
%Vss(p);
%pause();
%Wxx numerico
Wxx_num(1) = 0;
for i = 2:M(1/10)
    Wxx_num(i) = (Crank(i-1,tFin)-2*Crank(i,tFin)+Crank(i+1,tFin))*100;
end
Wxx_num(M(1/10)+1) = 2.8;
close all
plot(gridSpace,Wxx_num);