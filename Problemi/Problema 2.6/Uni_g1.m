% Funzione g1(x)
g1 = @(x) (1 + exp(-x.^2)) ./ (x.^2 + 3);

% Derivata calcolata manualmente (già semplificata)
g1_deriv = @(x) ((x.^2 + 3) .* (-2 .* x .* exp(-x.^2)) - (1 + exp(-x.^2)) .* (2 .* x)) ./ (x.^2 + 3).^2;

% Valori di x su [0,1]
x_vals = linspace(0, 1, 1000);

% Valori del modulo della derivata
g1_deriv_vals = abs(g1_deriv(x_vals));

% Controllo se sempre < 1
if all(g1_deriv_vals < 1)
    disp('La derivata è sempre minore di 1 su [0,1]');
else
    disp('Attenzione: la derivata NON è sempre < 1 su [0,1]');
end

% Valore massimo
fprintf('Valore massimo di |g_1''(x)| su [0,1]: %.6f\n', max(g1_deriv_vals));

% Grafico
figure;
plot(x_vals, g1_deriv_vals, 'b', 'LineWidth', 2); hold on;
yline(1, '--r', 'LineWidth', 2);
xlabel('x'); ylabel('|g''_1(x)|');
title('Modulo della derivata di g_1(x)');
legend('|g_1''(x)|', 'y = 1');
grid on;
