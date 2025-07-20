# Analisi del Metodo SOR per la Risoluzione di Sistemi Lineari

## Obiettivo

Analizzare il comportamento del metodo iterativo **SOR (Successive Over-Relaxation)** applicato alla risoluzione del sistema lineare:

\[
Ax = b
\]

e confrontarne le prestazioni con il metodo di **Gauss-Seidel**, determinando il valore ottimale del parametro di rilassamento \( \omega \) che minimizza il raggio spettrale \( \rho(G_\omega) \).

---

## Parte 1: Impostazione del Problema

Si considera il seguente sistema lineare:

\[
A =
\begin{bmatrix}
1 & -\frac{1}{4} & \frac{1}{3} & 0 \\
-1 & 2 & 0 & \frac{1}{2} \\
2 & 1 & 3 & -\frac{1}{3} \\
-1 & -2 & -4 & 7
\end{bmatrix}, \quad
b = \begin{bmatrix}
1 \\ 0 \\ -2 \\ 0
\end{bmatrix}
\]

Si parte da un vettore iniziale nullo \( x^{(0)} = 0 \).

---

## Parte 2: Decomposizione della Matrice

La matrice \( A \) viene decomposta come:

\[
A = D - E - F
\]

dove:
- \( D \) è la matrice diagonale contenente gli elementi diagonali di \( A \),
- \( E \) è la parte inferiore stretta di \( A \) (negata),
- \( F \) è la parte superiore stretta di \( A \) (negata).

Questa decomposizione è fondamentale per l'implementazione di metodi iterativi classici come Jacobi, Gauss-Seidel e SOR.

---

## Parte 3: Metodo SOR e Matrice di Iterazione

Nel metodo SOR, la formula iterativa è definita come:

\[
x^{(k+1)} = (D - \omega E)^{-1}[(1 - \omega)D + \omega F]x^{(k)} + \omega(D - \omega E)^{-1}b
\]

La matrice di iterazione associata è:

\[
G_\omega = (D - \omega E)^{-1}[(1 - \omega)D + \omega F]
\]

Nel codice, viene riscritta in forma equivalente numericamente stabile:

\[
M = \frac{1}{\omega}D - E, \quad N = \frac{1 - \omega}{\omega}D + F, \quad G_\omega = M^{-1}N
\]

Si calcola quindi \( \rho(G_\omega) \), cioè il raggio spettrale, per valori di \( \omega \in (0,2) \) in modo da trovare il minimo.

---

## Parte 4: Ottimizzazione di \( \omega \)

Tramite una discretizzazione fine di \( \omega \in (0,2) \), si determina:

\[
\omega_{\text{opt}} = \arg\min_{\omega \in (0,2)} \rho(G_\omega)
\]

Tale valore ottimizza la velocità di convergenza del metodo SOR.

---

## Parte 5: Confronto con Gauss-Seidel

Il metodo di **Gauss-Seidel** è un caso particolare del metodo SOR con \( \omega = 1 \). Si calcolano quindi i raggi spettrali:

- \( \rho(G_1) \): metodo Gauss-Seidel,
- \( \rho(G_{\omega_{\text{opt}}}) \): SOR ottimizzato.

---

## Parte 6: Iterazioni e Confronto

Vengono effettuate 10 iterazioni per entrambi i metodi:

- **Gauss-Seidel:**
\[
x^{(k+1)} = (D - E)^{-1}(F x^{(k)} + b)
\]

- **SOR con \( \omega_{\text{opt}} \):**
\[
x^{(k+1)} = M^{-1}(N x^{(k)} + b)
\]

---

## Parte 7: Valutazione dell’Errore

La soluzione esatta viene calcolata con:

\[
x = A^{-1}b
\]

e si confronta con le soluzioni approssimate ai vari step usando la norma euclidea:

\[
\|x^{(k)} - x\|_2
\]

Vengono infine plottati i grafici delle norme degli errori in scala logaritmica (semilogy) per evidenziare la **velocità di convergenza** dei due metodi.

---

## Conclusione

Il metodo SOR con \( \omega_{\text{opt}} \) converge significativamente più rapidamente rispetto al metodo di Gauss-Seidel, come confermato:

- dal confronto dei raggi spettrali \( \rho(G) \),
- dalla diminuzione dell’errore nei primi 10 step.

La scelta ottimale del parametro di rilassamento \( \omega \) è quindi cruciale per l’efficienza del metodo SOR.
