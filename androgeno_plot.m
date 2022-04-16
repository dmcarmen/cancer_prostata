function [] = androgeno_plot(dias_con, dias_sin, t_fin)
    arguments
        dias_con double = 180    %dias de tratamiento
        dias_sin double = 180    %dias sin tratamiento
        t_fin double = 1500
    end

    close all

    andro = @(t) A(t,dias_con,dias_sin);
    t = 0:1:t_fin;

    for i = t
        and(i+1) = andro(i);
    end
    plot(t,and)
    legend('A: andrógeno');
    title('Andrógeno ', ['con ' num2str(dias_con) ' dias con tratamiento y ' num2str(dias_sin) ' dias sin tratamiento'])
    xlabel('Tiempo (dias)')
end