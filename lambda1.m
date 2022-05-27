function [res] = lambda1(q)
    %LAMBDA1 Funcion lambda para el primer par de variables. 
    % Representa los ratios de transicion debido a mutaciones de las celulas de X1 (AD) a las de X2 (AI).
    % Params:
    %   q: valor de Q(t).
    % return:
    %   res: valor de la funcion.
    
    c = [0.00016 0.00012];
    K = [0.8 1.7];
    res = c(1) * (K(1)^3/(K(1)^3 + q^3));
end
