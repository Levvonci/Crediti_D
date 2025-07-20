function matlab()
    format long;
    % Risoluzione del Problema 2.6 completo
    
    % Definizione delle funzioni g1 e g2
    g1 = @(x) (1 + exp(-x.^2)) ./ (x.^2 + 3);
    g2 = @(x) cos(x);
    
    % Intervalli
    intervallo1 = [0, 1];
    intervallo2 = [0, pi/3];
    
    % Valori di epsilon da considerare
    epsilons = 10.^(-1:-1:-10);
    
    % Punto di innesco iniziale
    x0_1 = 0.5; % Puoi scegliere un valore qualsiasi in [0,1]
    x0_2 = pi/6; % Un valore qualsiasi in [0, pi/3]
    
    % Numero massimo di iterazioni
    Nmax = 1000;
    
    % Grafici
    figure;
    fplot(g1, intervallo1, 'r', 'LineWidth', 1.5);
    hold on;
    fplot(@(x) x, intervallo1, 'b--', 'LineWidth', 1.5);
    title('Grafico g1(x) e y=x');
    legend('g1(x)', 'y=x');
    grid on;
    
    figure;
    fplot(g2, intervallo2, 'r', 'LineWidth', 1.5);
    hold on;
    fplot(@(x) x, intervallo2, 'b--', 'LineWidth', 1.5);
    title('Grafico g2(x) e y=x');
    legend('g2(x)', 'y=x');
    grid on;
    
    % Tabelle risultati
    disp('--- Risultati per g1(x) ---');
    tabella1 = calcolaTabella(g1, x0_1, epsilons, Nmax);
    disp(tabella1);
    
    disp('--- Risultati per g2(x) ---');
    tabella2 = calcolaTabella(g2, x0_2, epsilons, Nmax);
    disp(tabella2);
end

function risultati = calcolaTabella(g, x0, epsilons, Nmax)
    % Calcolo approssimazioni alpha_e per diverse soglie epsilon
    
    n = length(epsilons);
    risultati = table('Size',[n 4],'VariableTypes',{'double','double','double','double'}, ...
                      'VariableNames',{'Epsilon','Alpha_e','NumIterazioni','ErroreModulo'});
    
    for i = 1:n
        epsilon = epsilons(i);
        [alpha_e, K, errore] = punto_fisso(g, x0, epsilon, Nmax);
        risultati.Epsilon(i) = epsilon;
        risultati.Alpha_e(i) = alpha_e;
        risultati.NumIterazioni(i) = K;
        risultati.ErroreModulo(i) = errore;
    end
end

function [xk, k, errore] = punto_fisso(g, x0, epsilon, Nmax)
    % Metodo di iterazione funzionale (punto fisso)
    
    xk = x0;
    for k = 1:Nmax
        xk1 = g(xk);
        if abs(xk1 - xk) <= epsilon
            xk = xk1;
            errore = abs(xk - g(xk));
            return
        end
        xk = xk1;
    end
    % Se non si soddisfa la condizione entro Nmax iterazioni
    errore = abs(xk - g(xk));
end
