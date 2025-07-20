function [Sn] = CavaSimp(a,b,f,n)
    % input:
    % a = estremo sinistro dell'intervallo
    % b = estremo destro dell'intervallo
    % f = funzione integrabile su [a,b]
    % n = numero naturale >=1
    %
    % output:
    % Sn = approssimazione dell'integrale su [a,b] della
    %      funzione f ottenuta mediante la formula di
    %      Cavalieri-Simpson di ordine n
    
    h=(b-a)/n;
    %s1 rappresenta la prima sommatoria della formula
    s1=0;
    for j=1:(n-1)
        xj=a+(j*h);
        s1=s1+f(xj);
    end

    %s2 rappresenta la seconda sommatoria della formula
    s2=0;
    for j=0:(n-1)
        x0=a+(j*h);
        x1=a+((j+1)*h);
        x2=(x0+x1)/2;
        s2=s2+f(x2);
    end

    %Sn rappresenta la formula di Cavalieri-Simpson di ordine n
    Sn=(h/6)*(f(a)+f(b)+2*s1+4*s2);
end