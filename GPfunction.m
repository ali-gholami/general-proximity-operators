function [phi,phi1,omega]=GPfunction(p,q,order)
% GPfunction generates the general potential function for a given p (0<p<=2)
% and q >= -1
% the outputs are phi and s which are returned as function handles
% Ali Gholami,
% University of Tehran, Iran.
if q==0
    phi =@(x) log(abs(x).^p+1);
else
    phi =@(x) 1/q*(1 - (abs(x).^p+1).^(-q));
end
phi1 =@(x) (p*x.^(p-1)./(x.^p+1).^(q+1)); % the first derivative
phi2 =@(x) ((p-1)*(x.^p+1) - p*(q+1)*x.^p).*(p*x.^(p-2))./(x.^p+1).^(q+2); % the second derivative

if strcmp(order, 'first')
    omega =@(x,tau) tau*phi1(x);
else 
    omega =@(x,tau) tau*phi1(x)./(1 + tau*phi2(x));
end
