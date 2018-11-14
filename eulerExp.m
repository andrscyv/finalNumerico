function [w,t] = eulerExp(f,y0,n,t0,tf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
h = (tf-t0)/n;
t  = zeros(n+1);
w = zeros(n+1);
w(1) = y0;
t(1) = t0;
for k=1:n
    t(k+1) = (k)*h+t0;
    w(k+1) = w(k)+h*f(t(k),w(k));
end

end

