# Problema 2.4 - Metodo SOR, Gauss-Seidel, Confronto di Convergenza

---

# Punto a - Tracciare \( \omega \mapsto \rho(G_\omega) \) e determinare \( \omega_{\text{opt}} \)

**Obiettivo:**  
Studiare come cambia il raggio spettrale \( \rho(G_\omega) \) al variare di \( \omega \), per trovare il valore ottimale che minimizza \( \rho(G_\omega) \).

**Procedura:**
- Costruiamo la matrice di iterazione:
  \[
  G_\omega = M^{-1}N
  \]
  dove
  \[
  M = \frac{1}{\omega} D - E, \quad N = \frac{1-\omega}{\omega} D + F
  \]
  e \( D \), \( E \), \( F \) sono le parti diagonale, inferiore e superiore di \( A \).
- Per ogni \( \omega \in (0,2) \) calcoliamo \( \rho(G_\omega) \), il massimo valore assoluto degli autovalori di \( G_\omega \).
- Individuiamo \( \omega_{\text{opt}} \) come il valore di \( \omega \) che minimizza \( \rho(G_\omega) \).
- Tracciamo il grafico \( \omega \mapsto \rho(G_\omega) \) per visualizzare il comportamento.

### Codice
```matlab
% Dati del problema
A = [1, -1/4, 1/3, 0;
    -1, 2, 0, 1/2;
     2, 1, 3, -1/3;
    -1, -2, -4, 7];
b = [1; 0; -2; 0];
m = 1000;
x0 = zeros(4,1);
epsilon = 1e-10;
Nmax = 10000;

% Decomposizione
D = diag(diag(A));
E = -tril(A, -1);
F = -triu(A, 1);

% Variazione di omega
omega_values = linspace(0.01, 2, 2*m);
rho_values = zeros(size(omega_values));

for i = 1:length(omega_values)
    omega = omega_values(i);
    M = (1/omega)*D - E;
    N = ((1-omega)/omega)*D + F;
    G_omega = M \ N;   % Matrice di iterazione
    rho_values(i) = max(abs(eig(G_omega))); % raggio spettrale
end

% Trova omega ottimale
[rho_min, idx_min] = min(rho_values);
omega_opt = omega_values(idx_min);

% Plot
figure;
plot(omega_values, rho_values, 'LineWidth', 2);
xlabel('\omega');
ylabel('\rho(G_\omega)');
title('Raggio spettrale in funzione di \omega');
grid on;
hold on;
plot(omega_opt, rho_min, 'ro', 'MarkerFaceColor','r');
legend('\rho(G_\omega)', 'omega_{opt}');
hold off;

fprintf('Omega ottimale: %.4f\n', omega_opt);
fprintf('Rho minimo: %.4f\n', rho_min);

```
### Grafico della funzione

![[assets/GraficoFunzione.png|center|500]]

### Risultato
**$\omega \ ottimale: 1.1240$**
**$\rho \ minimo: 0.2075$**

---

# Punto b - Calcolare \( \rho(G_\omega) \) per \( \omega = 1 \) e \( \omega = \omega_{\text{opt}} \)

**Obiettivo:**  
Confrontare il raggio spettrale:
- Quando \( \omega = 1 \) (corrispondente al metodo Gauss-Seidel).
- Quando \( \omega = \omega_{\text{opt}} \) (ottimale).

**Procedura:**
- Calcolare \( G_\omega \) e \( \rho(G_\omega) \) nei due casi.
- Verificare che \( \rho(G_{\omega_{\text{opt}}}) < \rho(G_1) \), ossia che il metodo ottimizzato converge più velocemente.

### Codice
```matlab
% Per omega = 1 (Gauss-Seidel)
omega1 = 1;
M1 = (1/omega1)*D - E;
N1 = ((1-omega1)/omega1)*D + F;
G1 = M1 \ N1;
rho1 = max(abs(eig(G1)));

fprintf('Rho(G) per omega = 1: %.4f\n', rho1);

% Per omega = omega_opt (già calcolato sopra)
Mopt = (1/omega_opt)*D - E;
Nopt = ((1-omega_opt)/omega_opt)*D + F;
Gopt = Mopt \ Nopt;
rho_opt = max(abs(eig(Gopt)));

fprintf('Rho(G) per omega_opt: %.4f\n', rho_opt);
```
### Risultato

**$\rho(G_\omega) \ per \ \omega= 1: 0.4604$**
**$\rho(G_\omega) \ per \ \omega_{opt} : 0.2075$**

---

# Punto c - Prime 10 iterazioni: confronto tra Gauss-Seidel e SOR con \( \omega_{\text{opt}} \)

**Obiettivo:**  
Osservare praticamente la velocità di convergenza dei due metodi in 10 iterazioni.

**Procedura:**
- Metodo Gauss-Seidel: risolvo
  \[
  x^{(k+1)} = (D - E)^{-1}(F x^{(k)} + b)
  \]
- Metodo SOR con \( \omega_{\text{opt}} \): risolvo
  \[
  x^{(k+1)} = M^{-1}(N x^{(k)} + b)
  \]
- Calcolare la norma \( \| x^{(k)} - x \|_2 \) ad ogni iterazione.
- Confrontare graficamente l'andamento dell'errore.

### Codice
```matlab
% Metodo Gauss-Seidel (omega = 1)
[x_GS, ~, ~] = metodo_SOR(A, b, 1, epsilon, x0, 10);

% Metodo SOR con omega_opt
[x_SORopt, ~, ~] = metodo_SOR(A, b, omega_opt, epsilon, x0, 10);

% Soluzione esatta
x_exact = A \ b;

% Iterazioni per confronto
X_GS = zeros(4,10);
X_SOR = zeros(4,10);

x_current = x0;
for k = 1:10
    x_current = (D - E) \ (F*x_current + b);
    X_GS(:,k) = x_current;
end

x_current = x0;
for k = 1:10
    M = (1/omega_opt)*D - E;
    N = ((1-omega_opt)/omega_opt)*D + F;
    x_current = M \ (N*x_current + b);
    X_SOR(:,k) = x_current;
end

% Tabelle dei risultati
disp('Iterazioni Gauss-Seidel:');
disp(X_GS);

disp('Iterazioni SOR con omega_opt:');
disp(X_SOR);

% Norme della differenza dalla soluzione esatta
norm_GS = vecnorm(X_GS - x_exact);
norm_SOR = vecnorm(X_SOR - x_exact);

figure;
semilogy(1:10, norm_GS, '-o', 1:10, norm_SOR, '-x');
xlabel('Numero di iterazioni');
ylabel('Norma errore ||x_k - x||_2');
legend('Gauss-Seidel', 'SOR \omega_{opt}');
title('Confronto convergenza tra Gauss-Seidel e SOR');
grid on;
```
### Tabelle

- **Iterazioni Gauss-Seidel**

| Iterazione | 1       | 2       | 3       | 4       | 5       | 6       | 7       | 8       | 9       | 10      |
|:-----------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|
| $x_1$         | 1.0000  | 1.6250  | 1.9495  | 2.0982  | 2.1667  | 2.1983  | 2.2128  | 2.2195  | 2.2225  | 2.2240  |
| $x_2$         | 0.5000  | 0.9554  | 1.1530  | 1.2443  | 1.2863  | 1.3056  | 1.3145  | 1.3186  | 1.3205  | 1.3214  |
| $x_3$         | -1.5000 | -2.1319 | -2.4299 | -2.5670 | -2.6301 | -2.6591 | -2.6725 | -2.6787 | -2.6815 | -2.6828 |
| $x_4$         | -0.5714 | -0.7132 | -0.7806 | -0.8116 | -0.8259 | -0.8324 | -0.8355 | -0.8369 | -0.8375 | -0.8378 |

- **Iterazioni SOR con $\omega_{opt}$**

| Iterazione | 1       | 2       | 3       | 4       | 5       | 6       | 7       | 8       | 9       | 10      |
|:-----------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|
| $x_1$         | 1.1240  | 1.8470  | 2.1444  | 2.2051  | 2.2213  | 2.2243  | 2.2250  | 2.2251  | 2.2252  | 2.2252  |
| $x_2$         | 0.6316  | 1.1819  | 1.2829  | 1.3150  | 1.3203  | 1.3218  | 1.3220  | 1.3221  | 1.3221  | 1.3221  |
| $x_3$         | -1.8282 | -2.4483 | -2.6329 | -2.6723 | -2.6816 | -2.6834 | -2.6838 | -2.6839 | -2.6839 | -2.6839 |
| $x_4$         | -0.7908 | -0.7983 | -0.8358 | -0.8363 | -0.8380 | -0.8380 | -0.8380 | -0.8380 | -0.8380 | -0.8380 |

### Confronto Convergenza

![[assets\ConfrConv.png|center|500]]

---

### Conclusioni

Dal confronto risulta evidente che:
- Il metodo SOR con paramentro ottimale $\omega_{opt}$ converge molto piú velocemente rispetto al metodo Gauss-Seidel con $\omega = 1$
- $\rho(G_{\omega_{opt}})$ é inferiore a $\rho(G_1)$
- Giá dopo poche iterazioni, il metodo SOR ottimizzato presenta un errore significativamente minore

---

# Codice Finale
```matlab
clc;
clear;
close all;

% Dati del problema
A = [1, -1/4, 1/3, 0;
    -1, 2, 0, 1/2;
     2, 1, 3, -1/3;
    -1, -2, -4, 7];
b = [1; 0; -2; 0];
x0 = zeros(4,1);
epsilon = 1e-10;
Nmax = 1000;
m = 1000; % Per la discretizzazione di omega

% Decomposizione matrice
D = diag(diag(A));
E = -tril(A, -1);
F = -triu(A, 1);

% (a) Calcolo rho(G_omega) per omega in (0,2)
omega_values = linspace(0.01, 2, 2*m);
rho_values = zeros(size(omega_values));

for i = 1:length(omega_values)
    omega = omega_values(i);
    M = (1/omega)*D - E;
    N = ((1-omega)/omega)*D + F;
    G_omega = M \ N;
    rho_values(i) = max(abs(eig(G_omega)));
end

% Trova omega ottimale
[rho_min, idx_min] = min(rho_values);
omega_opt = omega_values(idx_min);

% Grafico
figure;
plot(omega_values, rho_values, 'LineWidth', 2);
xlabel('\omega');
ylabel('\rho(G_\omega)');
title('Raggio spettrale in funzione di \omega');
grid on;
hold on;
plot(omega_opt, rho_min, 'ro', 'MarkerFaceColor','r');
legend('\rho(G_\omega)', 'omega_{opt}');
hold off;

fprintf('Omega ottimale: %.4f\n', omega_opt);
fprintf('Raggio spettrale minimo: %.4f\n', rho_min);

% (b) Calcolo rho(G) per omega=1 e omega=omega_opt
% Omega = 1 (Gauss-Seidel)
M1 = (1/1)*D - E;
N1 = ((1-1)/1)*D + F;
G1 = M1 \ N1;
rho1 = max(abs(eig(G1)));
fprintf('Rho(G) per omega=1: %.4f\n', rho1);

% Omega = omega_opt
Mopt = (1/omega_opt)*D - E;
Nopt = ((1-omega_opt)/omega_opt)*D + F;
Gopt = Mopt \ Nopt;
rho_opt = max(abs(eig(Gopt)));
fprintf('Rho(G) per omega_opt: %.4f\n', rho_opt);

% (c) 10 iterazioni Gauss-Seidel e SOR(omega_opt)

% Metodo Gauss-Seidel
X_GS = zeros(4,10);
x_current = x0;
for k = 1:10
    x_current = (D - E) \ (F*x_current + b);
    X_GS(:,k) = x_current;
end

% Metodo SOR con omega_opt
X_SOR = zeros(4,10);
x_current = x0;
for k = 1:10
    M = (1/omega_opt)*D - E;
    N = ((1-omega_opt)/omega_opt)*D + F;
    x_current = M \ (N*x_current + b);
    X_SOR(:,k) = x_current;
end

% Soluzione esatta
x_exact = A\b;

% Errori
norm_GS = vecnorm(X_GS - x_exact);
norm_SOR = vecnorm(X_SOR - x_exact);

% Tabelle
disp('Iterazioni Gauss-Seidel:');
disp(X_GS);

disp('Iterazioni SOR con omega_opt:');
disp(X_SOR);

% Grafico errori
figure;
semilogy(1:10, norm_GS, '-o', 'LineWidth', 2);
hold on;
semilogy(1:10, norm_SOR, '-x', 'LineWidth', 2);
xlabel('Numero di iterazioni');
ylabel('Norma errore ||x_k - x||_2');
legend('Gauss-Seidel', 'SOR \omega_{opt}');
title('Confronto convergenza tra Gauss-Seidel e SOR');
grid on;
hold off;
```