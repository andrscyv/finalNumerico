function [V, x] = FwdTCS_Calor(D, X, T, CF, CI, h, k)

a = X(1); b = X(2);
t0 = T(1); tf = T(2);

M = ceil( (b - a)/h ) + 1;
x = linspace(a, b, M);

N = ceil( (tf - t0)/k ) + 1;

% La solucion al tiempo inicial
V(:, 1) = CI(x)';

sigma = D*k/h^2;

aux = sigma*ones(1,M-3);
A = (1 - 2*sigma)*eye(M-2) + diag(aux,1) + diag(aux,-1);
F = zeros(M-2,1);

for i = 1:N
    F(1) = sigma*V(1,i); F(M-2) = sigma*V(M,i);
    V(2:M-1,i+1) = A*V(2:M-1,i) + F;

    % ATENCION: Esta condicion de frontera es tipo Dirichlet.
    t = t0 + i*k;
    V(1,i+1) = CF{1}(t);
    V(M,i+1) = CF{2}(t);
end

