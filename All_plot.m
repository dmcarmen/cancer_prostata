function [] = All_plot(dias_con, dias_sin, t_fin)
    %ALL_PLOT Grafica con las 5 EDOs.
    % Param:
    %   dias_con: dias con tratamiento de supresion androgeno. Por defecto
    %    60 dias.
    %   dias_sin: dias de descanso del tratamiento de supresion androgeno. 
    %    Por defecto 60 dias.
    %   t_fin: tiempo en dias de fin del plot. Por defecto 500.
    arguments
        dias_con double = 60    %dias de tratamiento
        dias_sin double = 60    %dias sin tratamiento
        t_fin double = 500
    end

    close all

    % Datos necesarios para el caculo del vector inicial. 
    alpha = 14.9/15;
    beta = 1;           %caso 1-4
    P0 = 15;
    t_ini = 0;
    %Q1(0) Q2(0) X1(0) X2(0) P(0)
    Var0 = [0.4 0.4 alpha*beta*P0 (1-alpha)*beta*P0 P0]; 

    S = @(t,y) dAlldt(t,y,dias_con, dias_sin);
    h=10e-3;

    % Grafica usando ode45.
    %ode45(S, [t_ini t_fin], Var0);
    [T, Y] = ode45(S, t_ini:h:t_fin, Var0);
    %[T, Y] = ode15s(S, [t_ini:h:t_fin], Var0);
    figure(1)
    plot(T,Y);
    legend('Q_1: cuota celular de AD', 'Q_2: cuota celular de AI', 'X_1: células AD', 'X_2: células AI', 'PSA');
    %lgd = legend({'Q_1: cuota celular de AD', 'Q_2: cuota celular de AI', 'X_1: células AD', 'X_2: células AI', 'PSA'}, 'Location','northwest');
    %lgd.NumColumns = 2;
    title('Solución utilizando ode45', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    % Grafica usando mi funcion de Euler explicito.
    %n = (t_fin-t_ini)*2;
    [~, y] = mi_euler(S, t_ini, t_fin, Var0, h);
    figure(2)
    plot(T,y(5,:)) % PSA
    legend('PSA')
    title('Solución utilizando función de Euler propia', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    % Grafica cone el error entre ambos metodos
    figure(3)
    plot(T,abs(Y-y')) % error entre ambos metodos
    legend('Error Q_1', 'Error Q_2', 'Error X_1', 'Error X_2', 'Error PSA')
    title('Error entre ode45 y función de Euler propia', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')

    %plot(t, y(i,:)); plot una sola variable en este orden Q1, Q2, X1, X2, P
    %plot(t, y(1:2,:)) ;plot Qis
    %plot(t, y(3:4,:)) ;plot Xis
    %plot(t, y(5,:)) ;plot P
    
end