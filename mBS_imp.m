function [W] = mBS_imp(Ix, It, M, N, bs)
% Proposito: esquema Euler implicito (Backward Time - Central Space)
% In : Ix ... intervalo del espacio [a, b]
% It ... intervalo del tiempo [0, T]
% M ... numero de pasos en el espacio
% N ... numero de pasos en el tiempo
% bs ... (estructura) caracteristicas de BS con
% bs.r ... taza de intereses
% bs.sigma ... coeficiente de volatibilidad
% bs.fc ... (function_handle) condicion terminal
% bs.bcL ... (function_handle) condicion iquierda
% bs.bcR ... (function_handle) condicion derecha
% Out: W ... (M+1)x(N+1) matriz con la aproximacion numerica
A = -(bs.r - ((bs.sigma)^2 )/ 2 )/ log(2);
B = - (((bs.sigma)^2) / 2)/ (log(2)^2);

h = (Ix(2) - Ix(1))/M;
k = (It(2) - It(1))/N;

lambda = k/h;
mu = k / h^2;

theta = B*mu - A*lambda/2;
phi = A*lambda/2 + B*mu;
gamm = 1+ bs.r*k - 2*B*mu;

D = zeros(M,M);
D = diag(theta*ones(M-1,1),-1) + diag(gamm*ones(M,1),0) + diag(phi*ones(M-1,1),1);
%D(1,1) = D(1,1)*bs.bcL(0);
D(M,M-2) = 1/(2*h);
D(M,M-1) = -2/h;
D(M,M) = 3/(2*h);

W = zeros(M+1,N+1);
%W(:,1) = bs.fc(Ix(1) + (0:M)*h);
for j =1:M+1
    W(j,1) = bs.fc(Ix(1)+(j-1)*h);
end

%plot(linspace(-2,2,M+1),W(:,1))

aux = zeros(M,1);
f_n = zeros(M,1);
%f_n(1,1) = theta*bs.bcL(It(1)+k);
%D = inv(D);
for i=2:N+1
    f_n(1,1) = theta*bs.bcL(It(1)+(i-1)*k);
    aux(1:M-1,1) = W(2:M,i-1);
    aux(M,1) = bs.bcR(It(1)+(i-1)*k);
    W(2:M+1,i) = D\(aux-f_n);
    %W(2:M+1,i) = D * (aux-f_n);
    %verifi = norm( D*W(2:M+1,i) -(aux-f_n))
    W(1,i) = bs.bcL(It(1)+(i-1)*k);
    
    
end 


%{
 for m =2:M
    for n = 2:N+1
       Err(m-1,n) =  W(m,n-1)-theta*W(m-1,n)-gamm*W(m,n)-phi*W(m+1,n);
    end
end
unos = ones(1,M-1);
for n = 1:N+1
   [val]  = min(abs(Err(:,n)))
   prom = unos*Err(:,n)/M-1
end 
%}


