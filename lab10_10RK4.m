close all; clear all;
f = @(t,y) 10*(1-y);
M = 70;
acum =0;
for k=1:M
tic
[w,t,h] = RK4(0,20,1/2,f,320);
acum =acum+ toc;
end
acum=acum/M
h
plot(t,w);