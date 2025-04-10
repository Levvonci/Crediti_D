function [xK, K, err] = punto_fisso(g, x0, epsilon, Nmax)
    % Metodo del punto fisso per risolvere x = g(x)
    %
    % Input:
    %   g       - funzione anonima (e.g., @(x) cos(x))
    %   x0      - punto iniziale
    %   epsilon - soglia di precisione
    %   Nmax    - numero massimo di iterazioni
    %
    % Output:
    %   xK  - valore dell'approssimazione finale
    %   K   - numero di iterazioni effettuate
    %   err - modulo dell'errore |x_K - g(x_K)|
    
    xK = x0;
    for K = 1:Nmax
        x_next = g(xK);
        if abs(x_next - xK) <= epsilon
            xK = x_next;
            err = abs(xK - g(xK));
            return;
        end
        xK = x_next;
    end
    % Se non si è soddisfatta la condizione di arresto
    err = abs(xK - g(xK));
end
