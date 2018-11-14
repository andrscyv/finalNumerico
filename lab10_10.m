close all; clear all; clc;
f = @(t,y) 10*(1-y);
M = 70;
acum =0;
for k=1:M
tic
[w,t] = eulerExp(f,1/2,320,0,20);
acum =acum+ toc;
end
acum=acum/M
plot(t,w);