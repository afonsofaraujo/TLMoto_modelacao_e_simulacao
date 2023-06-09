clear; clc; close all;
% parametros da simulação
dt = 0.001;         % intervalos
ts = 10;            % tempo da simulação
steps = ts./dt;
% constantes
R = 0.3;
m = 20;
g = 9.81;
k = 9800;
c = 1000;
I = [0.4875, 0, 0;
    0, 2.1000, 0;
    0, 0, 0.4875];
pressao_pneus = 2;
% Alocacao
model.x = zeros(3, steps);
model.v = zeros(3, steps);
model.a = zeros(3, steps);
model.euler = zeros(3, steps);
model.euler_dot = zeros(3, steps);
model.euler_dotdot = zeros(3, steps);
% initial conditions
model.x(:,1) = [0; 0; R];
model.v(:,1) = [30; 0; 0];
model.a(:,1) = [0; 0; 0];
model.euler(:,1) = [0; 0; 10];
model.euler_dot(:,1) = [0; 0; 30/R];
model.euler_dotdot(:,1) = [0; 0; 0];
% loop
for i = 1:1:steps
    % show progress
    if i == 0.25*steps
        disp("25%")
    elseif i == 0.5*steps
        disp("50%")
    elseif i == 0.75*steps
        disp("75%")
    end

    % limitar a velocidade
    velocidade = sqrt(model.v(1,i)^2 + model.v(2,i)^2 + model.v(3,i)^2);
    if  velocidade < 0.001
        disp("Velocidade atingiu 0 m/s - valor limite")
        break
    elseif velocidade > 10000
        disp("Velocidade atingiu 10000 m/s - valor limite")
        break
    end

    % limitar o camber
    if  abs(rad2deg(model.euler(2,i))) > 90
        disp("Camber atingiu 90° - valor limite")
        break
    end
    
    % definir d
    if velocidade <= 165/3.6
        f_w = 0.0085+0.018/pressao_pneus + velocidade^2*1.59e-6/pressao_pneus;
    else
        f_w = 0.018/pressao_pneus + velocidade^2*2.91e-6/pressao_pneus;
    end
    d = f_w*R;
    

    % momentos aplicados
    if (4 < i*dt) && (i*dt < 4.5)
        M_steer_val = 0.000001;
    else
        M_steer_val = 0;
    end
    if (2 < i*dt) && (i*dt < 2.5)
        M_acel_val = 50;
    else
        M_acel_val = 0;
    end


    % Matrizes de rotacao
    R_psi = [cos(model.euler(1,i)), sin(model.euler(1,i)), 0;
            -sin(model.euler(1,i)), cos(model.euler(1,i)), 0;
            0, 0, 1];
    R_phi = [1, 0, 0;
            0, cos(model.euler(2,i)), sin(model.euler(2,i));
            0, -sin(model.euler(2,i)), cos(model.euler(2,i))];
    R_theta = [cos(model.euler(3,i)), 0, sin(model.euler(3,i));
            0, 1, 0;
            -sin(model.euler(3,i)), 0, cos(model.euler(3,i))];
   
    % Variaveis
    if model.x(3,i) < R
        Normal = k*(R*cos(model.euler(2,i))-model.x(3,i)) + c*(-model.euler_dot(2,i)*R*sin(model.euler(2,i))-model.v(3,i));
    else
        Normal = 0;
    end
    slip = -(model.v(1,i)*cos(model.euler(1,i)) - model.v(2,i)*sin(model.euler(1,i)) - model.euler_dot(3,i)*R) / (model.v(1,i)*cos(model.euler(1,i)) - model.v(2,i)*sin(model.euler(1,i)));
    sideslip = atan((model.v(1,i)*sin(model.euler(1,i)) + model.v(2,i)*cos(model.euler(1,i))) / (model.v(1,i)*cos(model.euler(1,i)) - model.v(2,i)*sin(model.euler(1,i))));
    F_long = 20*slip*Normal;
    F_lat = Normal*(20*sideslip + 1.3*model.euler(2,i));


    % Forcas
    P = [0; 0; -m*g];
    N = [0; 0; Normal];
    F_long_vec = [F_long; 0; 0];
    F_lat_vec = [0; F_lat; 0];

    % Equilibrio de forcas - atualizacao da aceleracao
    model.a(:,i+1) = (1/m).*(P + N + R_psi'*F_lat_vec + R_psi'*F_long_vec);


    % Momentos
    r_f_lat = [0; 0; -R];                                
    M_f_lat = cross(r_f_lat, R_psi*R_phi*F_lat_vec);     % movel                  
    r_f_long = [0; 0; -R];                               
    M_f_long = cross(r_f_long, R_psi*R_phi*F_long_vec);  % movel
    r_f_normal = [d; 0; -R];
    M_f_normal = cross(r_f_normal, R_phi*N);             % movel

    M_steer = [M_steer_val; 0,; 0];                      % movel
    M_acel = [0; M_acel_val; 0];                         % movel
   

    % Equilibrio de momentos
    Aux_1 = [-cos(model.euler(2,i))*sin(model.euler(3,i)), cos(model.euler(3,i)), 0;
            sin(model.euler(2,i)), 0, 1;
            cos(model.euler(2,i))*cos(model.euler(3,i)), sin(model.euler(3,i)), 0];
    Aux_2 = [model.euler_dot(2,i)*sin(model.euler(2,i))*sin(model.euler(3,i))-model.euler_dot(3,i)*cos(model.euler(2,i))*cos(model.euler(3,i)), -model.euler_dot(3,i)*sin(model.euler(3,i)), 0;
            model.euler_dot(2,i)*cos(model.euler(2,i)), 0, 0;
            -model.euler_dot(2,i)*sin(model.euler(2,i))*cos(model.euler(3,i))-model.euler_dot(3,i)*cos(model.euler(2,i))*sin(model.euler(3,i)), model.euler_dot(3,i)*cos(model.euler(3,i)), 0];
    
    model.euler_dotdot(:,i+1) = (I*Aux_1)\(M_f_lat + M_f_long + M_f_normal + R_theta*M_steer + M_acel - I*Aux_2*model.euler_dot(:,i));
    
    
    % Atualizar CM e taxas de variacao
    model.v(:,i+1) = model.v(:,i) + dt*model.a(:,i);
    model.x(:,i+1) = model.x(:,i) + dt*model.v(:,i);
    model.euler_dot(:,i+1) = model.euler_dot(:,i) + dt*model.euler_dotdot(:,i);
    model.euler(:,i+1) = model.euler(:,i) + dt*model.euler_dot(:,i);
end
disp("Simulação terminada")


figure(1);
t = 1:steps;
subplot(3, 3, 1);
plot(t, model.x(1,t))
title("x")
subplot(3, 3, 2);
plot(t, model.x(2,t))
title("y")
subplot(3, 3, 3);
plot(t, model.x(3,t))
title("z")
subplot(3, 3, 4);
plot(t, model.v(1,t))
title("v_x")
subplot(3, 3, 5);
plot(t, model.v(2,t))
title("v_y")
subplot(3, 3, 6);
plot(t, model.v(3,t))
title("v_z")
subplot(3, 3, 7);
plot(t, model.a(1,t))
title("a_x")
subplot(3, 3, 8);
plot(t, model.a(2,t))
title("a_y")
subplot(3, 3, 9);
plot(t, model.a(3,t))
title("a_z")



figure(2);
% plot euler
t = 1:steps;
subplot(2, 3, 1);
plot(t, rad2deg(model.euler(1,t)),"Color",'k')
title("\psi")
ylabel("Graus [^\circ]")
subplot(2, 3, 2);
plot(t, rad2deg(model.euler(2,t)),"Color",'b')
title("\phi")
ylabel("Graus [^\circ]")
subplot(2, 3, 3);
plot(t, rad2deg(model.euler(3,t)),"Color",'r')
title("\theta")
ylabel("Graus [^\circ]")
% plot euler_dot
subplot(2, 3, 4);
plot(t, model.euler_dot(1,t),"Color",'k')
hold on
plot(t, model.euler_dot(2,t),"Color",'b')
plot(t, model.euler_dot(3,t),"Color",'r')
ylabel("[Rads/s]")
title("model.euler\_dot")
hold off
% plot euler_dotdot
subplot(2, 3, 6);
plot(t, model.euler_dotdot(1,t),"Color",'k')
hold on
plot(t, model.euler_dotdot(2,t),"Color",'b')
plot(t, model.euler_dotdot(3,t),"Color",'r')
ylabel("[Rads/s^2]")
title("model.euler\_dot\_dot")
hold off
