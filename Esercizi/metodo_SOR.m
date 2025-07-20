function [x, K, res_norm] = metodo_SOR(A, b, omega, epsilon, x0, Nmax)
    % Metodo SOR per risolvere il sistema lineare Ax = b
    % 
    % INPUT:
    % - A: matrice del sistema (n x n)
    % - b: vettore termine noto (n x 1)
    % - omega: parametro di rilassamento (omega ≠ 0)
    % - epsilon: soglia di precisione per il residuo
    % - x0: vettore di innesco iniziale (n x 1)
    % - Nmax: numero massimo di iterazioni
    %
    % OUTPUT:
    % - x: approssimazione della soluzione
    % - K: numero di iterazioni effettuate
    % - res_norm: norma 2 del residuo all'ultima iterazione
    
    %n = length(b);
    x = x0;
    
    % Decomposizione della matrice
    D = diag(diag(A));       % Parte diagonale
    E = -tril(A, -1);         % Parte inferiore (senza diagonale)
    F = -triu(A, 1);          % Parte superiore (senza diagonale)
    
    % Matrice del metodo
    M = (1/omega)*D - E;
    N = ((1-omega)/omega)*D + F;
    
    % Iterazione
    for K = 1:Nmax
        x_new = M \ (N*x + b);    % \ = risoluzione sistema lineare
        
        % Residuo
        r = b - A*x_new;
        res_norm = norm(r, 2);
        
        % Controllo condizione di arresto
        if res_norm <= epsilon
            x = x_new;
            return
        end
        
        % Aggiorna iterato
        x = x_new;
    end
end

