clc;
clear;
close all;

% Dati del problema
A = [1, -1/4, 1/3, 0;
    -1, 2, 0, 1/2;
     2, 1, 3, -1/3;
    -1, -2, -4, 7];
b = [1; 0; -2; 0];
x0 = zeros(4,1);
epsilon = 1e-10;
Nmax = 1000;
m = 1000; % Per la discretizzazione di omega

% Decomposizione matrice
D = diag(diag(A));
E = -tril(A, -1);
F = -triu(A, 1);

% (a) Calcolo rho(G_omega) per omega in (0,2)
omega_values = linspace(0.01, 2, 2*m);
rho_values = zeros(size(omega_values));

for i = 1:length(omega_values)
    omega = omega_values(i);
    M = (1/omega)*D - E;
    N = ((1-omega)/omega)*D + F;
    G_omega = M \ N;
    rho_values(i) = max(abs(eig(G_omega)));
end

% Trova omega ottimale
[rho_min, idx_min] = min(rho_values);
omega_opt = omega_values(idx_min);

% Grafico
figure;
plot(omega_values, rho_values, 'LineWidth', 2);
xlabel('\omega');
ylabel('\rho(G_\omega)');
title('Raggio spettrale in funzione di \omega');
grid on;
hold on;
plot(omega_opt, rho_min, 'ro', 'MarkerFaceColor','r');
legend('\rho(G_\omega)', 'omega_{opt}');
hold off;

fprintf('Omega ottimale: %.4f\n', omega_opt);
fprintf('Raggio spettrale minimo: %.4f\n', rho_min);

% (b) Calcolo rho(G) per omega=1 e omega=omega_opt
% Omega = 1 (Gauss-Seidel)
M1 = (1/1)*D - E;
N1 = ((1-1)/1)*D + F;
G1 = M1 \ N1;
rho1 = max(abs(eig(G1)));
fprintf('Rho(G) per omega=1: %.4f\n', rho1);

% Omega = omega_opt
Mopt = (1/omega_opt)*D - E;
Nopt = ((1-omega_opt)/omega_opt)*D + F;
Gopt = Mopt \ Nopt;
rho_opt = max(abs(eig(Gopt)));
fprintf('Rho(G) per omega_opt: %.4f\n', rho_opt);

% (c) 10 iterazioni Gauss-Seidel e SOR(omega_opt)

% Metodo Gauss-Seidel
X_GS = zeros(4,10);
x_current = x0;
for k = 1:10
    [x_current, ~, ~] = metodo_SOR(A, b, 1, epsilon, x_current, 1);
    X_GS(:,k) = x_current;
end

% Metodo SOR con omega_opt
X_SOR = zeros(4,10);
x_current = x0;
for k = 1:10
    [x_current, ~, ~] = metodo_SOR(A, b, omega_opt, epsilon, x_current, 1);
    X_SOR(:,k) = x_current;
end

% Soluzione esatta
x_exact = A\b;

% Errori
norm_GS = vecnorm(X_GS - x_exact);
norm_SOR = vecnorm(X_SOR - x_exact);

% Tabelle
disp('Iterazioni Gauss-Seidel:');
disp(X_GS);

disp('Iterazioni SOR con omega_opt:');
disp(X_SOR);

% === Calcolo norma infinito ===
norminf_GS = zeros(1, 10);
norminf_SOR = zeros(1, 10);

for k = 1:10
    norminf_GS(k) = norm(X_GS(:,k) - x_exact, inf);
    norminf_SOR(k) = norm(X_SOR(:,k) - x_exact, inf);
end

% Grafico errori
figure;
semilogy(1:10, norminf_GS, '-o', 'LineWidth', 2);
hold on;
semilogy(1:10, norminf_SOR, '-x', 'LineWidth', 2);
xlabel('Numero di iterazioni');
ylabel('Norma errore ||x_k - x||_\infty');
legend('Gauss-Seidel', 'SOR \omega_{opt}');
title('Confronto convergenza tra Gauss-Seidel e SOR');
grid on;
hold off;

% Stampa norme infinito
disp('Norma infinito Gauss-Seidel:');
disp(norminf_GS);

disp('Norma infinito SOR con omega_opt:');
disp(norminf_SOR);