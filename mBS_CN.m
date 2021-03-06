function [ W ] = mBS_CN ( Ix, It, M, N, bs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = -(bs.r - ((bs.sigma)^2 )/ 2 )/ log(2);
B = - (((bs.sigma)^2) / 2)/ (log(2)^2);

h = (Ix(2) - Ix(1))/M;
k = (It(2) - It(1))/N;

gamm = -1/k - bs.r/2 + B/h^2;
phi = -A/(4*h) - B/(2*h^2);
theta = A / (4*h) - B/(2*h^2);
delta = -1/k + bs.r/2 - B/h^2;

W = zeros(M+1,N+1);

BGrande = zeros(M);
BGrande = diag(theta*ones(M-1,1),-1) + diag(gamm*ones(M,1),0) + diag(phi*ones(M-1,1),1);
BGrande(M,M-2) = 1/(2*h);
BGrande(M,M-1) = -2/h;
BGrande(M,M) = 3/(2*h);

AGrande = zeros(M);
AGrande = diag(-theta*ones(M-1,1),-1) + diag(delta*ones(M,1),0) + diag(-phi*ones(M-1,1),1);

for j =1:M+1
    W(j,1) = bs.fc(Ix(1)+(j-1)*h);
end

f_n = zeros(M,1);
D = zeros(M,1);
BGrande = inv(BGrande);
for i = 2:N+1
    f_n(1,1) = -theta*(bs.bcL(It(1)+k*(i-2)) - bs.bcL(It(1)+k*(i-1)));
    D = AGrande*W(2:M+1,i-1);
    D(M,1) = bs.bcR(It(1)+(i-1)*k);
    W(2:M+1,i) = BGrande*(D+f_n);
    W(1,i) = bs.bcL(It(1)+(i-1)*k);
end

