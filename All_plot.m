function [] = All_plot(dias_con, dias_sin, t_fin)
    %ALL_PLOT Grafica con las 5 EDOs.
    % Param:
    %   dias_con: dias con tratamiento de supresion androgeno. Por defecto
    %    60 dias.
    %   dias_sin: dias de descanso del tratamiento de supresion androgeno. 
    %    Por defecto 60 dias.
    %   t_fin: tiempo en dias de fin del plot. Por defecto 500.
    arguments
        dias_con double = 190    %dias de tratamiento
        dias_sin double = 190    %dias sin tratamiento
        t_fin double = 1000
    end

    close all

    %% Datos necesarios para el caculo del vector inicial. 
    alpha = 14.9/15;
    beta = 1;           %caso 1-4
    P0 = 15;
    t_ini = 0;
    %Q1(0) Q2(0) X1(0) X2(0) P(0)
    Var0 = [0.4 0.4 alpha*beta*P0 (1-alpha)*beta*P0 P0]; 

    %% Ecuaciones diferenciale, division y calculo.
    S = @(t,y) dAlldt(t,y,dias_con, dias_sin);
    h=10e-3;

    [T, Y] = ode45(S, t_ini:h:t_fin, Var0);
    [~, y] = mi_euler(S, t_ini, t_fin, Var0, h);

    %% Graficas usando ode45.
    % Q1 y Q2.
    figure(1)
    plot(T,Y(:,1:2));
    legend('Q_1: cuota celular de AD', 'Q_2: cuota celular de AI');
    %lgd = legend({'Q_1: cuota celular de AD', 'Q_2: cuota celular de AI', 'X_1: celulas AD', 'X_2: celulas AI', 'PSA'}, 'Location','northwest');
    %lgd.NumColumns = 2;
    title('Solucion utilizando ode45', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    % X1, X2 y PSA.
    figure(2)
    plot(T,Y(:,3:5))
    legend('X_1: celulas AD', 'X_2: celulas AI', 'PSA');
    title('Solucion utilizando ode45', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    %% Graficas usando mi funcion de Euler explicito. 
    % Q_1 y Q_2.
    figure(3)
    plot(T,y(1:2,:)) 
    legend('Q_1: cuota celular de AD', 'Q_2: cuota celular de AI')
    title('Solucion utilizando funcion de Euler propia', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    
    % X_1, X_2 y PSA.
    figure(4)
    plot(T,y(3:5,:)) 
    legend('X_1: celulas AD', 'X_2: celulas AI', 'PSA')
    title('Solucion utilizando funcion de Euler propia', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    %% Graficas comparando resultados
    % Q1 y Q2.
    figure(5)
    plot(T,y(1:2,:));
    hold on
    plot(T,Y(:,1:2));
    legend('Q_1: cuota celular de AD ode45', 'Q_2: cuota celular de AI  ode45', ...
        'Q_1: cuota celular de AD Euler', 'Q_2: cuota celular de AI Euler');
    %lgd = legend({'Q_1: cuota celular de AD', 'Q_2: cuota celular de AI', 'X_1: celulas AD', 'X_2: celulas AI', 'PSA'}, 'Location','northwest');
    %lgd.NumColumns = 2;
    title('Solucion utilizando ambos metodos', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    hold off

    % X1, X2 y PSA.
    figure(6)
    plot(T, y(3:5,:))
    hold on
    plot(T,Y(:,3:5))
    legend('X_1: celulas AD ode45', 'X_2: celulas AI ode45', 'PSA ode45', ...
        'X_1: celulas AD Euler', 'X_2: celulas AI Euler', 'PSA Euler');
    title('Solucion utilizando ambos metodos', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    hold off

    %% Grafica con el error entre ambos metodos
    figure(7)
    semilogy(T,abs(Y-y')) % error absoluto entre ambos metodos
    legend('Error Q_1', 'Error Q_2', 'Error X_1', 'Error X_2', 'Error PSA')
    title('Semilogy Error entre ode45 y funcion de Euler propia', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    %% Grafica del androgeno
    figure(8)
    for i = 1:length(T)
        and(i) = A(T(i), dias_con, dias_sin);
    end
    plot(T,and)
    legend('A: androgeno');
    title('Androgeno ', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
end
