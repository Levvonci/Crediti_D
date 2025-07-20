% Definizione della funzione f(x)
f = @(x) 1 ./ (x .* log(x));

% Intervallo di integrazione
a = 2;
b = 5;

% Calcolare il valore esatto dell'integrale I
I_exact = integral(f, a, b);

% Vettore dei valori di n
n_values = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560];

% Inizializzare vettori per gli errori
error_trapezi = zeros(size(n_values));
error_simpson = zeros(size(n_values));

% Creiamo la tabella
fprintf('%-8s %-15s %-15s %-15s %-15s\n', 'n', 'I_n', 'S_n', '|I_n - I|', '|S_n - I|');
fprintf('%s\n', repmat('-', 1, 68));

for i = 1:length(n_values)
    n = n_values(i);
    
    % Approssimazione tramite la formula dei trapezi
    I_n = Trapezi(a, b, n, f);
    
    % Approssimazione tramite la formula di Cavalieri-Simpson
    S_n = CavaSimp(a, b, f, n);
    
    % Calcolare gli errori assoluti
    error_trapezi(i) = abs(I_n - I_exact);
    error_simpson(i) = abs(S_n - I_exact);
    
    % Visualizzare i risultati
    fprintf('%-8d %-15.6e %-15.6e %-15.6e %-15.6e\n', n, I_n, S_n, error_trapezi(i), error_simpson(i));
end