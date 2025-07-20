function [xK, K, err] = fixed_point_iteration(g, epsilon, x0, Nmax)
% fixed_point_iteration risolve l'equazione x = g(x) usando il metodo di iterazione funzionale
% Input:
%   g       - funzione handle per g(x)
%   epsilon - soglia di precisione per il criterio di arresto
%   x0      - punto di innesco
%   Nmax    - numero massimo di iterazioni
% Output:
%   xK      - ultima approssimazione calcolata
%   K       - numero di iterazioni effettuate
%   err     - errore finale |xK - g(xK)|

x_prev = x0; % primo valore iniziale

for K = 1:Nmax
    xK = g(x_prev); % calcola il nuovo valore
    if abs(xK - x_prev) <= epsilon
        err = abs(xK - g(xK));
        return
    end
    x_prev = xK; % aggiorna per l'iterazione successiva
end

% Se raggiungo Nmax iterazioni senza soddisfare la condizione
err = abs(xK - g(xK));
end
