function [] = All_plot(dias_con, dias_sin, t_fin)
    %ALL_PLOT Grafica con las 5 EDOs.
    % Param:
    %   dias_con: dias con tratamiento de supresion androgeno. Por defecto 200 dias.
    %   dias_sin: dias de descanso del tratamiento de supresion androgeno. Por defecto 200 dias.
    %   t_fin: tiempo en dias de fin del plot. Por defecto 1800.
    arguments
        dias_con double = 200    %dias de tratamiento
        dias_sin double = 200    %dias sin tratamiento
        t_fin double = 1800
    end

    close all  

    %% Datos necesarios para el caculo del vector inicial. 
    alpha = 14.9/15;
    beta = 1;
    
    P0 = 15;
    %Var0 = [Q1(0) Q2(0) X1(0) X2(0) P(0)]
    Var0 = [0.4 0.4 alpha*beta*P0 (1-alpha)*beta*P0 P0]; 

    %% Ecuaciones diferenciales, calculo y error.
    t_ini = 0;

    S = @(t,y) dAlldt(t,y,dias_con, dias_sin);
    h=1e-2;

    [T, Y] = ode45(S, t_ini:h:t_fin, Var0);
    [~, y] = mi_euler(S, t_ini, t_fin, Var0, h);

    abs_error = abs(Y-y');

    %% Graficas usando ode45.
    % Q1 y Q2.
    figure(1)
    plot(T,Y(:,1:2))

    axis([0 t_fin -inf inf])
    axis 'auto y'

    legend({'Q_1: cuota celular de AD', 'Q_2: cuota celular de AI'}, 'Location','best');
    title('Solucion utilizando ode45', [num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    fontsize(gca,scale=1.25)
    saveas(1,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Qs_ode45.png'])

    % X1, X2 y PSA.
    figure(2)
    plot(T,Y(:,3:5))

    axis([0 t_fin 0 20])
    axis 'auto y'

    legend({'X_1: celulas AD', 'X_2: celulas AI', 'PSA'}, 'Location','north');
    title('Solucion utilizando ode45', [ num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    fontsize(gca,scale=1.25)
    saveas(2,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Xs_PSA_ode45.png'])

    %% Graficas usando mi funcion de Euler explicito. 
    % Q_1 y Q_2.
    figure(3)
    plot(T,y(1:2,:))

    legend('Q_1: cuota celular de AD', 'Q_2: cuota celular de AI')
    title('Solucion utilizando el metodo de Euler', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    saveas(3,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Qs_Euler.png'])
    
    % X_1, X_2 y PSA.
    figure(4)
    plot(T,y(3:5,:))

    legend('X_1: celulas AD', 'X_2: celulas AI', 'PSA')
    title('Solucion utilizando el metodo de Euler', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    saveas(4,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Xs_PSA_Euler.png'])

    %% Graficas con el error entre ambos metodos
    % Q1 y Q2
    figure(5)
    semilogy(T,abs_error(:,1:2))

    axis([0 t_fin -inf inf])
    axis 'auto y'

    legend({'Error Q_1', 'Error Q_2'}, 'Location','best');
    title('Semilogy Error de Q_1 y Q_2 entre ode45 y metodo de Euler', [ num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    fontsize(gca,scale=1.1)
    saveas(5,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Qs_error.png'])

    % X1, X2 y PSA
    figure(6)
    semilogy(T,abs_error(:,3:5))

    axis([0 t_fin -inf inf])
    axis 'auto y'

    legend({'Error X_1', 'Error X_2', 'Error PSA'}, 'Location','best')
    title('Semilogy Error de X_1, X_2 y PSA entre ode45 y metodo de Euler', [ num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    fontsize(gca,scale=1.1)
    saveas(6,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_Xs_PSA_error.png'])

    %% Grafica del androgeno
    figure(7)
    for i = 1:length(T)
        and(i) = A(T(i), dias_con, dias_sin);
    end
    plot(T,and)

    axis([0 t_fin -inf inf])
    axis 'auto y'

    legend({'A: androgeno'}, 'Location','best');
    title('Androgeno ', [ num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
    fontsize(gca,scale=1.1)
    saveas(7,['imgs/' num2str(dias_con) '_' num2str(dias_sin) '_' num2str(t_fin) '_androgeno.png'])
end
