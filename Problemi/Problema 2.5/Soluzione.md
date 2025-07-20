# Caso 1 
$$f(x) = x^3 + 3x - 1 - e^{-x^2}, [a,b]=[0,1]$$

### Punto (a): Verifica che $f(a)f(b) < 0$

1. Calcoliamo $f(a)$ e $f(b)$:
    - $f(0) = 0^3 + 3(0) - 1 - e^{-0^2} = -1 - 1 = -2$,
    - $f(1) = 1^3 + 3(1) - 1 - e^{-1^2} = 1 + 3 - 1 - e^{-1} = 3 - e^{-1} \approx 2.63$.
2. Poiché $f(0) \cdot f(1) < 0$, $(\text{ risulta }-2\cdot2,63=-5,26)$ possiamo procedere.

### Punto (b): Grafico di $f(x)$ e verifica di uno zero unico

Tracciamo il grafico di $f(x)$ su $[0,1]$ con MATLAB per osservare che $f(x)$ ha un unico zero nell'intervallo $(0,1)$.

```matlab
f = @(x) x.^3 + 3.*x - 1 - exp(-x.^2);
x = linspace(0, 1, 1000); % 1000 punti nell'intervallo [0, 1]
plot(x, f(x), 'b-', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('f(x)');
title('Grafico di f(x) = x^3 + 3x - 1 - e^{-x^2}');
```

**Analisi:** Osservando il grafico, si nota che $f(x)$ è continuo e cambia segno una sola volta tra $0$ e $1$.

![[assets/Caso1.png|center|500]]

### Punto (c): Dimostrazione analitica che f(x) ha un unico zero

Usiamo il teorema di Bolzano e la monotonicità derivata dall'analisi di $f'(x)$:$$f'(x) = 3x^2 + 3 + 2x e^{-x^2}.$$
1. $f'(x) > 0$ per ogni $x \in [0, 1]$ (la funzione è strettamente crescente su $[0,1]$).
2. Poiché $f(x)$ è crescente e cambia segno in $[0,1]$, per il teorema di Bolzano esiste un unico zero $\zeta \in (0, 1)$.

### Punto (d): Tabella per $\varepsilon \in \{10^{-1}, 10^{-2}, \dots, 10^{-10}\}$

Abbiamo usato il **metodo di bisezione** per calcolare:

- L'approssimazione $\xi_\varepsilon$,
- Il numero di iterazioni $K_\varepsilon$,
- Il valore $f(\xi_\varepsilon)$.

La tabella dei risultati è la seguente:

| $\epsilon$           | $x_i$               | $K$  | $f(x_i)$                            |
| -------------------- | ------------------- | ---- | ----------------------------------- |
| $1.0 \cdot 10^{-1}$  | $0.531250000000000$ | $4$  | $-1.041995243049776 \cdot 10^{-2}$  |
| $1.0 \cdot 10^{-2}$  | $0.535156250000000$ | $7$  | $7.765312582933004 \cdot 10^{-3}$   |
| $1.0 \cdot 10^{-3}$  | $0.533691406250000$ | $10$ | $9.389559548024229 \cdot 10^{-4}$   |
| $1.0 \cdot 10^{-4}$  | $0.533477783203125$ | $14$ | $-5.586409047664276 \cdot 10^{-5}$  |
| $1.0 \cdot 10^{-5}$  | $0.533489227294922$ | $17$ | $-2.574612559369527 \cdot 10^{-6}$  |
| $1.0 \cdot 10^{-6}$  | $0.533489704132080$ | $20$ | $-3.542067064099541 \cdot 10^{-7}$  |
| $1.0 \cdot 10^{-7}$  | $0.533489793539047$ | $24$ | $6.211948844203619 \cdot 10^{-8}$   |
| $1.0 \cdot 10^{-8}$  | $0.533489782363176$ | $27$ | $1.007871253122516 \cdot 10^{-8}$   |
| $1.0 \cdot 10^{-9}$  | $0.533489780034870$ | $30$ | $-7.631157927789900 \cdot 10^{-10}$ |
| $1.0 \cdot 10^{-10}$ | $0.533489780180389$ | $34$ | $-8.550160579545718 \cdot 10^{-11}$ |


**Codice MATLAB:**

```matlab
a = 0; b = 1;
f = @(x) x.^3 + 3.*x - 1 - exp(-x.^2);
epsilon_values = 10.^(-1:-1:-10); % Tolleranze
results = zeros(length(epsilon_values), 3); % Preallocazione: [xi_eps, K_eps, f(xi_eps)]

for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    [xi, K, fx] = bisezione(a, b, f, epsilon);
    results(i, :) = [xi, K, fx];
end

% Mostra la tabella
disp('Tabella dei risultati:');
disp('epsilon        xi_eps         K_eps        f(xi_eps)');
disp(results);
```

---

# Caso 2 
$$f(x) = \cos x - x, [a,b]=[0,\pi]$$

### Punto (a): Verifica che $f(a) f(b) < 0$

1. Calcoliamo $f(a)$ e $f(b)$:
    - $f(0) = \cos(0) - 0 = 1$,
    - $f(\pi) = \cos(\pi) - \pi = -1 - \pi < 0$.
2. Poiché $f(0) \cdot f(\pi) < 0$, possiamo procedere.

### Punto (b): Grafico di $f(x)$ e verifica di uno zero unico

**Codice MATLAB:**

```matlab
f = @(x) cos(x) - x;
x = linspace(0, pi, 1000); % 1000 punti nell'intervallo [0, pi]
plot(x, f(x), 'r-', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('f(x)');
title('Grafico di f(x) = cos(x) - x');
```

**Analisi:** Il grafico mostra che $f(x)$ è continuo e cambia segno una sola volta tra $0$ e $\pi$.

![[assets/Caso2.png|center|500]]


### Punto (c): Dimostrazione analitica che $f(x)$ ha un unico zero

Usiamo $f'(x) = -\sin(x) - 1$:

1. $f'(x) < 0$ per ogni $x \in [0, \pi]$ (la funzione è strettamente decrescente su $[0, \pi]$).
2. Poiché $f(x)$ è decrescente e cambia segno in $[0, \pi]$, per il teorema di Bolzano esiste un unico zero $\zeta \in (0, \pi)$.

### Punto (d): Tabella per 
$$\varepsilon \in \{10^{-1}, 10^{-2}, \dots, 10^{-10}\}$$
Abbiamo usato il **metodo di bisezione** per calcolare:

- L'approssimazione $\xi_\varepsilon$,
- Il numero di iterazioni $K_\varepsilon$,
- Il valore $f(\xi_\varepsilon)$.

La tabella dei risultati è la seguente:

| $\epsilon$           | $x_i$               | $K$  | $f(x_i)$                            |
| -------------------- | ------------------- | ---- | ----------------------------------- |
| $1.0 \cdot 10^{-1}$  | $0.736310778185108$ | $5$  | $4.640347169851511 \cdot 10^{-3}$   |
| $1.0 \cdot 10^{-2}$  | $0.739378739760879$ | $9$  | $-4.914153002637534 \cdot 10^{-4}$  |
| $1.0 \cdot 10^{-3}$  | $0.738995244563908$ | $12$ | $1.504357420498703 \cdot 10^{-4}$   |
| $1.0 \cdot 10^{-4}$  | $0.739043181463529$ | $15$ | $7.021030579146270 \cdot 10^{-5}$   |
| $1.0 \cdot 10^{-5}$  | $0.739088122306924$ | $19$ | $-5.002583233437718 \cdot 10^{-6}$  |
| $1.0 \cdot 10^{-6}$  | $0.739085500757726$ | $22$ | $-6.151237084139893 \cdot 10^{-7}$  |
| $1.0 \cdot 10^{-7}$  | $0.739085173064076$ | $25$ | $-6.669162500028136 \cdot 10^{-8}$  |
| $1.0 \cdot 10^{-8}$  | $0.739085135028206$ | $29$ | $-3.034334783436066 \cdot 10^{-9}$  |
| $1.0 \cdot 10^{-9}$  | $0.739085133199558$ | $32$ | $2.611200144997383 \cdot 10^{-11}$  |
| $1.0 \cdot 10^{-10}$ | $0.739085133245275$ | $35$ | $-5.039924033667376 \cdot 10^{-11}$ |

**Codice MATLAB:**

```matlab
a = 0; b = pi;
f = @(x) cos(x) - x;
epsilon_values = 10.^(-1:-1:-10); % Tolleranze
results = zeros(length(epsilon_values), 3); % Preallocazione: [xi_eps, K_eps, f(xi_eps)]

for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    [xi, K, fx] = bisezione(a, b, f, epsilon);
    results(i, :) = [xi, K, fx];
end

% Mostra la tabella
disp('Tabella dei risultati:');
disp('epsilon        xi_eps         K_eps        f(xi_eps)');
disp(results);
```

---

# Codice Finale
```matlab
% Funzioni e intervalli definiti dal problema
f1 = @(x) x.^3 + 3*x - 1 - exp(-x.^2); % Prima funzione
a1 = 0; b1 = 1; % Intervallo [a, b] per f1

f2 = @(x) cos(x) - x; % Seconda funzione
a2 = 0; b2 = pi; % Intervallo [a, b] per f2

% Lista di epsilon
epsilons = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];

% Risoluzione per il primo caso
solve_case(f1, a1, b1, epsilons, 'f1(x) = x^3 + 3x - 1 - e^{-x^2}');

% Risoluzione per il secondo caso
solve_case(f2, a2, b2, epsilons, 'f2(x) = cos(x) - x');

% Funzione per risolvere ogni caso
function solve_case(f, a, b, epsilons, case_name)
    fprintf('\nSoluzione per %s:\n', case_name);
    
    % (a) Verifica che f(a)*f(b) < 0
    fa = f(a);
    fb = f(b);
    fprintf('(a) f(a)*f(b) = %.3f (segno opposto: %s)\n', fa * fb, ...
        string(fa * fb < 0));
    if fa * fb >= 0
        error('f(a) e f(b) devono avere segni opposti');
    end
    
    % (b) Tracciamento del grafico
    fprintf('(b) Tracciamento del grafico di f(x) su [%f, %f]\n', a, b);
    fplot(f, [a b]);
    hold on;
    grid on;
    plot(a, f(a), 'ro', 'DisplayName', 'f(a)');
    plot(b, f(b), 'bo', 'DisplayName', 'f(b)');
    xlabel('x'); ylabel('f(x)');
    title(['Grafico di f(x) - Caso ', case_name]);
    legend show;

    % (c) Dimostrazione analitica: fatta in modo separato (se necessario)

    % (d) Tabella dei risultati per vari epsilon
    fprintf('(d) Calcolo del metodo di bisezione per diverse tolleranze epsilon:\n');
    fprintf('epsilon       xi               K       f(xi)\n');
    fprintf('-------------------------------------------------\n');
    for epsilon = epsilons
        [xi, K, fx] = bisezione(a, b, f, epsilon);
        fprintf('%e   %.15f   %d   %.15e\n', epsilon, xi, K, fx);
    end
end
```
### Spiegazione Codice Finale
1. **Caso 1 e Caso 2**:
    - Si calcolano $f(a)$ e $f(b)$ per verificare che il prodotto è negativo.
    - Si tracciano i grafici per osservare il comportamento di $f(x)$.
    - Si riportano i risultati delle tabelle usando il metodo di bisezione.
2. **Funzione `bisezione`**:
    - Implementa il metodo di bisezione per trovare l'approssimazione di uno zero di una funzione continua su un intervallo $[a,b]$.
    - Restituisce l'approssimazione $\xi_\varepsilon$, il numero di iterazioni $K_\varepsilon$, e il valore $f(\xi_\varepsilon)$.
3. **Tabelle dei Risultati**:
    - Si stampano le tabelle per ogni caso, con i valori di $\epsilon$, $\xi_\varepsilon$,$K_\varepsilon$, e $f(\xi_\varepsilon)$.