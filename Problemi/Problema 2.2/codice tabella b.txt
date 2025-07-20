% Definizione della funzione f(x) = exp(x)
f = @(x) exp(x);

% Valore esatto dell'integrale I
I_exact = exp(1) - 1;

% Valori di epsilon
epsilon_values = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];

% Preallocazione dei risultati
n_values = zeros(size(epsilon_values));
I_n_values = zeros(size(epsilon_values));
errors = zeros(size(epsilon_values));

% Calcolo di n(epsilon), I_n e dell'errore
for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    n_values(i) = ceil(sqrt(exp(1)/(12*epsilon))); % n(epsilon) arrotondato verso l'alto
    n = n_values(i);
    
    % Calcolo di I_n con la funzione Trapezi
    I_n = Trapezi(0, 1, n, f); % Chiamata alla funzione Trapezi
    I_n_values(i) = I_n;
    
    % Calcolo dell'errore
    errors(i) = abs(I_exact - I_n);
end

% Tabella dei risultati
fprintf('%-10s %-12s %-20s %-15s\n', 'Epsilon', 'n(epsilon)', 'I_n', 'Errore');
for i = 1:length(epsilon_values)
    fprintf('%-10.1e %-12d %-20.10f %-15.10e\n', epsilon_values(i), n_values(i), I_n_values(i), errors(i));
end

function [app] = Trapezi(a, b, n, f)
    % Formula dei trapezi
    r = 0;
    h = (b - a) / n;
    for j = 1:(n - 1)
        r = r + f(a + j * h);
    end
    app = (((f(a) + f(b)) / 2) + r) * h;
end
