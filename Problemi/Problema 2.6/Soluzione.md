# Punto a - Tracciare i grafici di $y = g(x)$ e $y = x$

```matlab
% Dati del problema
g1 = @(x) (1 + exp(-x.^2)) ./ (x.^2 + 3);
g2 = @(x) cos(x);

% Intervalli
a1 = 0; b1 = 1;
a2 = 0; b2 = pi/3;

% Grafico g1(x) e y = x
figure;
fplot(g1, [a1 b1], 'LineWidth', 2); hold on;
fplot(@(x) x, [a1 b1], '--r', 'LineWidth', 2);
title('Grafico di g_1(x) e y=x');
legend('g_1(x)', 'y=x');
grid on;

% Grafico g2(x) e y = x
figure;
fplot(g2, [a2 b2], 'LineWidth', 2); hold on;
fplot(@(x) x, [a2 b2], '--r', 'LineWidth', 2);
title('Grafico di g_2(x) e y=x');
legend('g_2(x)', 'y=x');
grid on;
```

Si noti che dai grafici generati vedremo che $x = g(x)$ ha un'unica soluzione di $\alpha$ in entrambi i casi, ed essa corrisponde all'intersezione delle curve $y = g(x)$ e $y = x$



### Grafico $g_1(x)$
![[assets/GraficoFunzione1.png|center|500]]



### Grafico $g_2(x)$
![[assets/GraficoFunzione2.png|center|500]]

---
# Punto b - Dimostrazione analitica dell'unicitá

Dobbiamo verificare che:

- Le funzioni \( g_1(x) \) e \( g_2(x) \) sono di classe \( C^1 \) (cioè derivabili e con derivata continua);
- \( |g_1'(x)| < 1 \) e \( |g_2'(x)| < 1 \) nei rispettivi intervalli.

Se queste due condizioni sono soddisfatte, allora possiamo applicare il **Teorema del Punto Fisso**: l'equazione \( x = g(x) \) ha **un'unica soluzione** \( \alpha \) nell'intervallo dato.



### Calcolo di \( g_1'(x) \)

La funzione è definita da:

\[
g_1(x) = \frac{1 + e^{-x^2}}{x^2 + 3}
\]

Applicando la regola del quoziente, otteniamo:

\[
g_1'(x) = \frac{(x^2+3)(-2xe^{-x^2}) - (1+e^{-x^2})(2x)}{(x^2+3)^2}
\]

Questa derivata è continua su \([0,1]\) e possiamo verificare numericamente (ad esempio con MATLAB) che:

\[
|g_1'(x)| < 1 \quad \text{per ogni} \quad x \in [0,1]
\]



### Calcolo di \( g_2'(x) \)

La funzione è:

\[
g_2(x) = \cos(x)
\]

e derivando otteniamo:

\[
g_2'(x) = -\sin(x)
\]

Poiché \( |\sin(x)| \leq \sin\left(\frac{\pi}{3}\right) \approx 0.866 \) su \([0, \pi/3]\), abbiamo:

\[
|g_2'(x)| < 1 \quad \text{per ogni} \quad x \in [0, \pi/3]
\]



## Conclusione

In entrambi i casi:

- \( g(x) \) è continua e derivabile,
- \( |g'(x)| < 1 \) nell'intervallo specificato,

quindi, per il **Teorema del Punto Fisso di Banach**, l'equazione \( x = g(x) \) ha **esattamente una soluzione** \( \alpha \) in ciascun intervallo.

\[ \boxed{\text{Soluzione unica garantita!}} \]



### In parole semplici:

Per dimostrare l'unicità, abbiamo controllato che la derivata di \( g(x) \) è:

- Continua (\( g \in C^1 \)),
- Stretta in valore assoluto sotto 1 (\( |g'(x)| < 1 \)),

questo implica che \( g \) è una **contrazione**. Una funzione contrattiva ha sempre **un solo punto fisso**, che è la soluzione dell'equazione \( x = g(x) \).

---

# Punto c - Costruzione della tabella

Si crei una tabella per i valori di $\epsilon \in \{10^{-1},10^{-2}, \dots, 10^{-10}\}$.
Preso un valore $x_{0}$ casuale all'interno dell'intervallo (Esempio: 0.5 per entrambi i casi)

```matlab
% Funzioni
g1 = @(x) (1 + exp(-x.^2)) ./ (x.^2 + 3);
g2 = @(x) cos(x);

% Intervalli
a1 = 0; b1 = 1;
a2 = 0; b2 = pi/3;

% Punto di innesco
x0_1 = 0.5;
x0_2 = 0.5;

% Valori di epsilon
epsilons = 10.^(-(1:10));
Nmax = 1000; % per sicurezza

% Tabelle vuote
results1 = [];
results2 = [];

for i = 1:length(epsilons)
    eps = epsilons(i);
    
    % Primo problema (g1)
    [alpha1, K1, err1] = fixed_point_iteration(g1, eps, x0_1, Nmax);
    results1 = [results1; eps, x0_1, K1, abs(alpha1 - g1(alpha1))];
    
    % Secondo problema (g2)
    [alpha2, K2, err2] = fixed_point_iteration(g2, eps, x0_2, Nmax);
    results2 = [results2; eps, x0_2, K2, abs(alpha2 - g2(alpha2))];
end

% Visualizzazione tabella
disp('Tabella per g1(x)');
disp('epsilon      x0       K       |alpha - g(alpha)|');
disp(results1);

disp('Tabella per g2(x)');
disp('epsilon      x0       K       |alpha - g(alpha)|');
disp(results2);
```

### Risultati per $g_1(x)$
| $\epsilon$ | $\alpha_e$ | $NumIterazioni$ | $ErroreModulo$ |
| ---------- | ---------- | --------------- | -------------- |
| $0.1$      | $0.54732$  | $1$             | $0.019635$     |
| $0.01$     | $0.53591$  | $3$             | $0.0034297$    |
| $0.001$    | $0.53331$  | $6$             | $0.00024999$   |
| $0.0001$   | $0.5335$   | $9$             | $1.8216e-05$   |
| $1e-05$    | $0.53349$  | $11$            | $3.1778e-06$   |
| $1e-06$    | $0.53349$  | $14$            | $2.3156e-07$   |
| $1e-07$    | $0.53349$  | $16$            | $4.0397e-08$   |
| $1e-08$    | $0.53349$  | $19$            | $2.9437e-09$   |
| $1e-09$    | $0.53349$  | $22$            | $2.145e-10$    |
| $1e-10$    | $0.53349$  | $24$            | $3.7421e-11$   |



### Risultati per $g_2(x)$
| $\epsilon$ | $\alpha_e$ | $NumIterazioni$ | $ErroreModulo$ |
| ---------- | ---------- | --------------- | -------------- |
| $0.1$      | $0.69859$  | $4$             | $0.067167$     |
| $0.01$     | $0.73535$  | $10$            | $0.0062432$    |
| $0.001$    | $0.73874$  | $16$            | $0.00058305$   |
| $0.0001$   | $0.73905$  | $22$            | $5.4469e-05$   |
| $1e-05$    | $0.73908$  | $28$            | $5.0887e-06$   |
| $1e-06$    | $0.73908$  | $34$            | $4.7541e-07$   |
| $1e-07$    | $0.73909$  | $39$            | $6.5935e-08$   |
| $1e-08$    | $0.73909$  | $45$            | $6.1599e-09$   |
| $1e-09$    | $0.73909$  | $51$            | $5.7549e-10$   |
| $1e-10$    | $0.73909$  | $57$            | $5.3764e-11$   |

---

# Codice Finale
```matlab
function problema()
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
```

### Spiegazione Codice Finale
- **problema()** é il programma principale: definisce le funzioni, traccia i grafici, chiama il calcolo tabellare.
- **calcolaTabella()**  cicla sui valori di $\epsilon$ e salva tutti i risultati in una tabella MATLAB.
- **puntoFisso()** implementa esattamente il metodo di iterazione funzionale.
- Se il metodo converge in meno di **Nmax** iterazioni, restituisce il valore trovato; se no, dà comunque l'ultima iterazione.