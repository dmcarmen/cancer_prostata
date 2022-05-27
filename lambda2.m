function [res] = lambda2(q)
    %LAMBDA2 Funcion lambda para el segundo par de variables. 
    % Representa los ratios de transicion debido a mutaciones de las celulas de X2 (AI) a las de X1 (AD).
    % Params:
    %   q: valor de Q(t).
    % return:
    %   res: valor de la funcion.
    
    c = [0.00016 0.00012];
    K = [0.8 1.7];
    res = c(2) * (q^3/(K(2)^3 + q^3));
end
