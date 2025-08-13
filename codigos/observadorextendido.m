clc; clear; close all;
% Parámetros
g = 9.81; L = 10.0; u = 0; A = 1.0;   % A: amplitud de sin(3t)
% Condiciones iniciales
x1_0 = 0; x2_0 = 0;
% Simulación
h = 1e-3; t_final = 10;
N = round(t_final/h);
t = (0:N-1)'*h;
% Asignar memoria
x1 = zeros(N,1); x2 = zeros(N,1);
x1(1) = x1_0;   x2(1) = x2_0;
for k = 1:N-1
    u = sin(t(k)); 
    pint = -(g/L) * sin(x1(k));   % perturbación “interna”: -g/L sin(theta)
    pext =  0.1*sin(3*t(k));      % perturbación externa: sin(3t)
    dx1 = x2(k);
    dx2 = u + pint + pext;        % ω̇ = u − (g/L)sinθ + sin(3t)
    % Forward Euler
    x1(k+1) = x1(k) + h*dx1;
    x2(k+1) = x2(k) + h*dx2;
end

% ----- Observador extendido (ESO) -----
lam1 = 1; lam2 = 40; lam3 = 30;
x1h = zeros(N,1); x2h = zeros(N,1); zh = zeros(N,1);
for k = 1:N-1
    % --- ESO ---
    e = x1(k) - x1h(k);                % error de salida (y - y_hat)
    dx1h = x2h(k) + lam1*e;
    dx2h = u + zh(k) + lam2*e;
    dzh  = lam3*e;
    x1h(k+1) = x1h(k) + h*dx1h;
    x2h(k+1) = x2h(k) + h*dx2h;
    zh(k+1)  = zh(k)  + h*dzh;
end

plot(t, x1); hold on; plot(t, x2);

plot(t,x1h,'--')
plot(t,x2h,'--')

legend('\theta (rad)','\omega (rad/s)'); grid on; xlabel('t (s)');
title('Péndulo con entrada u y perturbaciones interna/externa');
