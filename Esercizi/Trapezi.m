function [app] = Trapezi(a,b,n,f)
    % input:
    % a = estremo sinistro dell'intervallo
    % b = estremo destro dell'intervallo
    % n = numero naturale >=1
    % f = funzione integrabile su [a,b]
    %
    % output:
    % app = approssimazione dell'integrale su [a,b] della 
    %    funzione f ottenuta mediante la formula dei trapezi 
    %    di ordine n
    r=0;
    h=(b-a)/n;
    for j=1:(n-1)
        r=r+f(a+j*h);
    end
    app=((f(a)+f(b))/2 + r)*h;
end