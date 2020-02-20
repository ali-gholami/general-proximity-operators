clc, clear; close all
% a demo for deconvolution problem 
%    argmin_x ||y - Gx||_2^2 + tau \sum \phi_q^p(x)
%   Ali Gholami,
% University of Tehran, Iran.
% The reference paper is:
%  A General Framework for Sparsity-Based Denoising and Inversion
% IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 59, NO. 11, NOVEMBER 2011
% by Ali Gholami and S. Mohammad Hosseini

addpath('./functions')
n = 2^8;
k = 5;
r = zeros(n,1);
r([20 50 130 150 230])= [1 -.7 .9 -1 .5];
w = ricker(10,0.002); % Ricker wavelet
sigma = 0.05;
randn('state',11)


N = length(r) + length(w) - 1;
G =@(x,N) real(ifft(fft(w,N).*fft(x,N))); % 
Gt =@(x,N) real(ifft(conj(fft(w,N)).*fft(x,N)));

y = G(r,N) + sigma*randn(N,1); % trace

%  regularization parameter
subplot(2,2,1)
plot(1:n, r,'k','LineWidth',2)
axis([1 n -1 1]);
title('original signal')
% xlabel('sample number')
set(gca,'fontsize',10)
subplot(2,2,2)
plot(1:N, y,'k','LineWidth',2)
axis([1 N -1 1]);
title('original signal')
% xlabel('sample number')
set(gca,'fontsize',10)
L= max(abs(fft(w)))^2;
p=1;
q=5;
tau=5;
[phi,dphi,omega] = GPfunction(p,q,'first'); % poentioal function
ye =@(tau) find_ye(phi,dphi,tau);         % jumping point
prox =@(y,tau) T_fun(y,omega,tau,ye(tau));% Proximity operator

M = N;
y = padarray(y,[N,0]);
N = length(y);
x = Gt(y,N);

t_old = 1;
x_old = x;
u = x_old;

for k=1:200
    gamma = k/(k*L*2+1);
    x = prox(u + gamma*Gt(y - G(u,N),N),tau*gamma);
    % FISTA update
    t = (1 + sqrt(1+4*t_old^2))/2;
    u = x + (t_old - 1)/t*(x - x_old);
    t_old  = t;
    x_old = x;
end
%     figure(1)
subplot(2,2,3)
plot(1:n, x(M+(1:n)),'k','LineWidth',2)
axis([1 n -1 1]);