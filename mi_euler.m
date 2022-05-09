function  [t, y] = mi_euler(f,t0,T,y0,h)
    % MI_EULER Resuelve un sistema de EDOs mediante Euler explicito.
    % Params:
    %   f: Y'=f(t,Y) matriz con la expresion de las EDOs.
    %   t0: tiempo inicial.
    %   T: tiempo final.
    %   y0: Y(t0)=y0, vector con el valor inicial. 
    %   h: tamanio de las particionnes de en el intervalo [t0,T].
    % return:
    %   t: los ts donde se ha evaluado la funcion 
    %      (es decir, [t0,T] dividido en N pasos).
    %   y: la matriz solucion.

    N=(T-t0)/h;     % Numero de divisiones
    t= t0:h:T;      % Vector de tiempos, tamanio N+1

    % Inicializamos la matriz a 0, excepto el valor incial
    y=zeros(length(y0),N+1); 
    y(:,1)=y0;

    for n=1:N
        % Euler explicito: y_{n+1}=y_n+h*f(t_n,y_n)
        y(:,n+1)=y(:,n)+h*f(t(n),y(:,n));
    end
end