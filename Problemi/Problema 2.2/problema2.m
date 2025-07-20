% Definizione della funzione e dell'intervallo
f = @(x) exp(x);
a = 0; b = 1;

% Calcolo delle approssimazioni con la formula dei trapezi
I2 = Trapezi(a, b, 2, f);    % n = 2
I4 = Trapezi(a, b, 4, f);    % n = 4
I8 = Trapezi(a, b, 8, f);    % n = 8
I16 = Trapezi(a, b, 16, f);  % n = 16

% Nodi e valori per l'interpolazione
H = [1/2, 1/4, 1/8, 1/16].^2;  % Quadrati dei passi
I = [I2, I4, I8, I16];          % Approssimazioni corrispondenti

% Valutazione del polinomio interpolante con ValPol
T = 0; % Valutiamo il polinomio in x = 0
P0 = ValPol(H, I, T);

% Valore esatto dell'integrale
I_exact = exp(1) - 1;

% Stampa dei risultati
fprintf('Valore esatto I: %.10f\n', I_exact);
fprintf('Valori calcolati con i trapezi:\n');
fprintf('I2 = %.10f, I4 = %.10f, I8 = %.10f, I16 = %.10f\n', I2, I4, I8, I16);
fprintf('Valore di p(0): %.10f\n', P0);

% Confronto tra I2, I4, I8, I16, p(0) e il valore esatto I
fprintf('\nConfronto:\n');
fprintf('Errore |I2 - I| = %.10e\n', abs(I2 - I_exact));
fprintf('Errore |I4 - I| = %.10e\n', abs(I4 - I_exact));
fprintf('Errore |I8 - I| = %.10e\n', abs(I8 - I_exact));
fprintf('Errore |I16 - I| = %.10e\n', abs(I16 - I_exact));
fprintf('Errore |p(0) - I| = %.10e\n', abs(P0 - I_exact));
