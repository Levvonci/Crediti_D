clear; clc; close all;
%--------------------------------------------------------------------------
% 1) DEFINIZIONE DEI NODI DI INTERPOLAZIONE E DELLA FUNZIONE DA INTERPOLARE
%--------------------------------------------------------------------------
% I nodi sono i valori: 0, 1/64, 4/64, 9/64, 16/64, 25/64, 36/64, 49/64, 64/64
x_nodes = [0, 1/64, 4/64, 9/64, 16/64, 25/64, 36/64, 49/64, 64/64]; % [cite: 2]
% Funzione di cui vogliamo fare l'interpolazione
y_nodes = sqrt(x_nodes); % [cite: 3]

%--------------------------------------------------------------------------
% 2) COSTRUZIONE DEL POLINOMIO DI INTERPOLAZIONE p(x)
%-------------------------------------------------------------------------- 
% La valutazione del polinomio avverrà direttamente tramite ValPol.

%--------------------------------------------------------------------------
% 3) CALCOLO DEL VETTORE DELLE DIFFERENZE p(x_i) - sqrt(x_i) PER 21 PUNTI
%--------------------------------------------------------------------------
% Creiamo 21 punti equispaziati tra 0 e 1
N = 21; % [cite: 5]
x_eval = linspace(0, 1, N);    % suddivisione di [0,1] in 21 punti [cite: 6]

% Valutiamo il polinomio interpolante nei 21 punti usando ValPol
% È necessario che ValPol.m sia accessibile (stessa cartella o nel path)
p_eval = ValPol(x_nodes, y_nodes, x_eval);

% Calcoliamo la funzione sqrt(x) negli stessi 21 punti
f_eval = sqrt(x_eval); % [cite: 7]
% Vettore delle differenze: p(x_i) - sqrt(x_i)
diff_vector = p_eval - f_eval; % [cite: 8]

% Visualizziamo il vettore delle differenze in colonna
disp('Vettore p(x_i) - sqrt(x_i) per i = 1,...,21 (calcolato con ValPol):'); % [cite: 9]
% Stampiamo i valori in tre blocchi separati
disp('Columns 1 through 9:'); % [cite: 10]
disp(diff_vector(1:9)); % [cite: 10]

disp('Columns 10 through 18:'); % [cite: 10]
disp(diff_vector(10:18)); % [cite: 10]

disp('Columns 19 through 21:'); % [cite: 11]
disp(diff_vector(19:21)); % [cite: 11]


%--------------------------------------------------------------------------
% 4) GRAFICI: CONFRONTO TRA sqrt(x) E p(x)
%--------------------------------------------------------------------------
figure;
% Tracciamo sqrt(x)
fplot(@(x) sqrt(x), [0, 1], 'LineWidth', 2); % [cite: 12]
hold on; grid on;
% Tracciamo p(x) usando ValPol
% ValPol necessita dei nodi (x_nodes, y_nodes) e dei punti di valutazione (x)
fplot(@(x_dynamic) ValPol(x_nodes, y_nodes, x_dynamic), [0, 1], 'LineWidth', 2);

% Aggiungiamo titolo, legenda e assi
title('Confronto tra $\sqrt{x}$ e il polinomio di interpolazione p(x) (con ValPol)', ...
      'Interpreter', 'latex');
xlabel('x'); ylabel('y'); % [cite: 14]
legend('$\sqrt{x}$', 'p(x) (ValPol)', 'Location','best', 'Interpreter', 'latex'); % [cite: 14]