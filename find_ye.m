function [ye xz] = find_ye(phi,dphi,tau)
% FIND_YE find the ye corresponding to jumpping point for a given potential
% function phi(x) which is supplied as a function handle. dphi(x) is
% another function handel which provide the derivative of phi(x). tau is a
% positive scaler (regularization parameter.
% find_ye returns the maximum root of x^2*r(x)
% written by Ali Gholami,
% University of Tehran, Iran.
r =@(x) 1 + 2*tau*(x.*dphi(x)-phi(x))./x.^2;
x1 = 1e-8*tau;
x2 = 1e3*tau;
xm = fminbnd(r,x1,x2);
if (r(xm) >= 1)
    xz = 0;
    ye = 0;
    return
elseif (r(xm) >= 0)
    xz = (1 - r(xm))*xm;
elseif (r(xm)< 0)
    xz = fzero(r,[xm x2]);
end
ye = xz + tau*dphi(xz);
