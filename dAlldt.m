function dydt = dAlldt(t,y, dias_con, dias_sin)
    %dAlldt Devuelve el sistema de EDOs. 
    % Params:
    %   t: escalar t del que dependen las funciones.
    %   y: vector de las propias funciones (Q1, Q2, X1, X2 y P).
    %   dias_con: dias con tratamiento de supresion androgeno.
    %   dias_sin: dias de descanso del tratamiento de supresion androgeno. 
    % return:
    %   dydt: vector con las ecuaciones dQ1/dt, dQ2/dt, dX1/dt, dX2/dt y dP/dt.
    
    % Parametros para dQi/dt
    vm = 0.27;
    qm = 5;
    vh = 4;
    b = 0.09;
    qmin = [0.363 0.153]; 
    mum = 0.027;

    % Parametros para dXi/dt
    %qmin = [0.363 0.153];
    %mum = 0.027;
    d = [0.016 0.017];

    % Parametros para dP/dt
    sigma = [0.02 0.28 0.34];
    delta = 0.08;
    rho = [1.3 1.1];
    
    % Inicializamos el vector de ecuaciones diferenciales
    dydt = zeros(5,1);
    
    %dydt(1) = dQ1/dt, y(1) = Q1.
    s1 = vm*((qm - y(1))/(qm - qmin(1)));
    s2 = A(t, dias_con, dias_sin)/(A(t, dias_con, dias_sin)+vh);
    s3 = mum*(y(1)-qmin(1));
    s4 = b*y(1);
    dydt(1) = s1*s2 - s3 - s4;
    
    %dydt(2) = dQ2/dt, y(2) = Q2.
    s1 = vm*((qm - y(2))/(qm - qmin(2)));
    %s2 = A(t, dias_con, dias_sin)/(A(t, dias_con, dias_sin)+vh);
    s3 = mum*(y(2)-qmin(2));
    s4 = b*y(2);
    dydt(2) = s1*s2 - s3 - s4;
    

    %dydt(3) = dX1/dt, y(3) = X1.
    s1 = mum*(1-(qmin(1)/y(1)))*y(3);
    s2 = d(1) * y(3);
    s3 = lambda1(y(1))*y(3);
    s4 = lambda2(y(2))*y(4);
    dydt(3) = s1 - s2 - s3 + s4;
    
    %dydt(4) = dX2/dt, y(4) = X2.
    s1 = mum*(1-(qmin(2)/y(2)))*y(4);
    s2 = d(2) * y(4);
    s3 = lambda2(y(2))*y(4);
    s4 = lambda1(y(1))*y(3);
    dydt(4) = s1 - s2 - s3 + s4;

    %dydt(5) = dP/dt, y(5) = P.
    % sigma(1) = sigma0, sigma(2) = sigma1, sigma(3) = sigma2 
    s1 = sigma(1)*(y(3) + y(4));
    s2 = y(3) * ((sigma(2)*y(1)^3)/(y(1)^3+rho(1)^3));
    s3 = y(4) * ((sigma(3)*y(2)^3)/(y(2)^3+rho(2)^3));
    s4 = delta*y(5);
    dydt(5) = s1 + s2 + s3 - s4;
end
