function [V] = ValPol(X, Y, T)
    % input:
    % X = vettore a componenti reali tutti distinti 
    % Y = vettore a componenti reali della stessa lunghezza di X
    % T = vettore contenente i punti in cui verrà calcolato p(x)
    % 
    % output:
    % V = vettore che contiene le valutazioni nei punti del vettore T del
    %     polinomio p(x) interpolante i valori Y sui nodi X.
    
    
    C =CNewton(X,Y);
    %C=vettore dei coefficienti di Newton
    V = zeros(1, length(T));
    
    for k=1:length(T)
        V(k)=RH(T(k),C,X);
    end
    %V=vettore dei polinomi calcolati con t

    function [f]=CNewton(X,Y)
        Z =zeros(length(Y));
         %Z=tabelle delle differenze divise
        for i=1:length(Y)
            Z(i,1)=Y(i);
        end
       % f(1)=Z(1,1);
        for j=2:length(Y)
            for i=j:length(Y)
                Z(i,j)=(Z(i,j-1)-Z(j-1,j-1))/(X(i)-X(j-1));
             %   if i==j
              %      f(j)=Z(i,j);
               % end
            end
        end
        f = diag(Z);
    end


    function [g]=RH(t,C,X)
        g = 0;
        n = length(C);
        for i = n:-1:1
            g = g * (t - X(i)) + C(i);
        end
    end
end
