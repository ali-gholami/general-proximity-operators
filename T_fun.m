function x = T_fun(y,omega,tau,ye)
ay = abs(y);
if ye
    x = (y - sign(y).*omega(ay,tau)).*(ay > ye);
    x(y==0) = 0;
else
    x = (ay - omega(ay,tau));
    x = (abs(x) + x)/2;
    x = sign(y).*x;
    x(y==0) = 0;
end
    
