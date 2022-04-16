function  [val]=A(t, dias_con, dias_sin)
    %ANDROGENO funcion que modela el cambio de androgeno. 
    % Comienza con una fase de tratamiento, seguida de una fase de descanso
    % indefinidamente. Modelada por la solucion de la EDO en Ideta model 
    % (supplementary).
    % Params:
    %   t: tiempo en el que calcular el valor de A(t)
    %   dias_con: dias con tratamiento de supresion androgeno. Opcional, 
    %    por defecto 60 dias.
    %   dias_sin: dias de descanso del tratamiento de supresion androgeno. 
    %    Opcional, por defecto 60 dias.
    %   a0: normal androgen concentration. Por defecto 20.
    %   gamma: androgen clearance rate. Por defecto 0.08
    % return:
    %   val: valor de A(t).
    arguments
        t double
        dias_con double = 60    %dias de tratamiento
        dias_sin double = 60    %dias sin tratamiento
        %a0 double = 20          %normal androgen concentration
        %gamma double = 0.08     %androgen clearance rate
    end

    a0 = 15;
    gamma = 0.08;

    ini0 = a0-0.5;  %constante inicial
    tf1 = 0;    %tiempo inicial

    % Bucle hasta encontrar si t esta en un periodo de tratamiento o
    % descanso y calcular el valor de A teniendo eso en cuenta.
    while true
        tf2 = tf1+dias_con; %tiempo al terminar tratamiento
        tf3 = tf2+dias_sin; %tiempo al terminar descanso
        b = ini0*exp(-gamma*(tf2-tf1))+1; %valor de A al terminar tratamiento
        
        if (t >= tf1) && (t < tf2)
            % de tf1 a tf2 con tratamiento
            val = ini0*exp(-gamma*(t-tf1))+0.5;
            break
        end
        if (t >=tf2) && (t < tf3)
            % de tf2 a tf3 sin tratamiento
            val = (b-a0)*exp(-gamma*(t-tf2))+a0;
            break
        end
        
        ini0 = b*exp(-gamma*(tf3-tf1))+a0-0.5;  %valor de A al terminar descanso
        tf1 = tf3;  %nuevo tiempo inicial de nueva ronda tratamiento-descanso
    end
end
