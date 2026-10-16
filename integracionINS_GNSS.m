% Integracion INS/GNSS Loosely Coupled
% BRANDO FRAIZ FREIRE
%% Setup inicial
clear; clc; close all;
fprintf('Empezando la simulacion...\n');

%% Parametros de la simulacion
T_obs = 300;        % tiempo total [s]
h = 0.01;           % paso de tiempo pequeño (para el INS)
T_s = 1.0;          % paso de tiempo grande (para el GNSS)

%% 1. Generar la trayectoria real (ground truth)
fprintf('Generando trayectoria...\n');

% Parametros para que la trayectoria no sea una linea recta
params_N.V_media = 10.0; params_N.sigma_v = 3.0; params_N.T_corr = 20.0;
params_E.V_media = 5.0;  params_E.sigma_v = 2.0; params_E.T_corr = 15.0;
params_U.V_media = 1.0;  params_U.sigma_v = 1.0; params_U.T_corr = 10.0;
params_N.sigma_ruido = sqrt(2 * params_N.sigma_v^2 / (params_N.T_corr));
params_E.sigma_ruido = sqrt(2 * params_E.sigma_v^2 / (params_E.T_corr));
params_U.sigma_ruido = sqrt(2 * params_U.sigma_v^2 / (params_U.T_corr));

num_puntos = T_obs / h + 1;
t_real = (0:h:T_obs)';
pos_enu_real = zeros(num_puntos, 3);
vel_enu_real = zeros(num_puntos, 3);

% Integrando para crear el camino
pos_E = 0; pos_N = 0; pos_U = 0;
v_fluct_E = 0; v_fluct_N = 0; v_fluct_U = 0;
for k = 1:num_puntos
    V_total_E = params_E.V_media + v_fluct_E;
    V_total_N = params_N.V_media + v_fluct_N;
    V_total_U = params_U.V_media + v_fluct_U;
    vel_enu_real(k, :) = [V_total_E, V_total_N, V_total_U];
    pos_enu_real(k, :) = [pos_E, pos_N, pos_U];
    pos_E = pos_E + h * V_total_E; pos_N = pos_N + h * V_total_N; pos_U = pos_U + h * V_total_U;
    v_fluct_E = v_fluct_E + h * (-v_fluct_E / params_E.T_corr + randn() * params_E.sigma_ruido/sqrt(h));
    v_fluct_N = v_fluct_N + h * (-v_fluct_N / params_N.T_corr + randn() * params_N.sigma_ruido/sqrt(h));
    v_fluct_U = v_fluct_U + h * (-v_fluct_U / params_U.T_corr + randn() * params_U.sigma_ruido/sqrt(h));
end

% Aceleracion real, la necesitamos para simular el INS
accel_enu_real = [zeros(1,3); diff(vel_enu_real)/h];

%% 2. Simular GNSS (con ruido)
fprintf('Simulando datos del GNSS...\n');

sigma_gnss_pos = 5.0;  % metros
sigma_gnss_vel = 0.5;  % m/s

% Muestreamos la trayectoria real y le metemos ruido
idx_muestreo = 1:round(T_s/h):num_puntos;
t_filtro = t_real(idx_muestreo);
num_pasos_filtro = length(t_filtro);

pos_gnss_enu = pos_enu_real(idx_muestreo, :) + randn(num_pasos_filtro, 3) * sigma_gnss_pos;
vel_gnss_enu = vel_enu_real(idx_muestreo, :) + randn(num_pasos_filtro, 3) * sigma_gnss_vel;

%% 3. Simular INS (con deriva)
fprintf('Simulando datos del INS...\n');

% Errores del acelerometro
sigma_acc_noise = 0.05;         % ruido normal del sensor
initial_bias_acc = [0.03, -0.02, 0.05]; % bias que tiene de fabrica
sigma_acc_drift = 0.001;        % Si es muy bajo, el filtro no corrige bien la deriva.

% Variables para el INS
pos_ins_only_enu = zeros(num_puntos, 3);
vel_ins_only_enu = zeros(num_puntos, 3);
bias_acc_real = zeros(num_puntos, 3);

% Empezamos bien, para ver como se degrada
pos_ins_only_enu(1,:) = pos_enu_real(1,:);
vel_ins_only_enu(1,:) = vel_enu_real(1,:);
bias_acc_real(1,:) = initial_bias_acc;

% Bucle de integracion del INS (Dead Reckoning)
for k = 2:num_puntos
    % El bias va cambiando solo con el tiempo (random walk)
    bias_acc_real(k,:) = bias_acc_real(k-1,:) + randn(1,3) * sigma_acc_drift * sqrt(h);
    
    % La medida del acelerometro es la real + bias + ruido
    accel_medida = accel_enu_real(k,:) + bias_acc_real(k,:) + randn(1,3) * sigma_acc_noise;
    
    % Integramos para sacar velocidad y posicion
    vel_ins_only_enu(k,:) = vel_ins_only_enu(k-1,:) + accel_medida * h;
    pos_ins_only_enu(k,:) = pos_ins_only_enu(k-1,:) + vel_ins_only_enu(k-1,:)*h + 0.5*accel_medida*h^2;
end

% Hacemos una copia de la solucion INS para ir corrigiendola
pos_ins_prop = pos_ins_only_enu;
vel_ins_prop = vel_ins_only_enu;

%% 4. Filtro de Kalman Loosely Coupled
fprintf('Ejecutando el filtro de Kalman...\n');

% Vector de estado del error (9 estados)
% x = [delta_p, delta_v, bias_acc]'
x_lc = zeros(9, 1); 

% Covarianza del estado P
P_lc = diag([10^2, 10^2, 10^2, ...       % incertidumbre inicial en pos
             1^2, 1^2, 1^2, ...         % incertidumbre inicial en vel
             (0.1)^2, (0.1)^2, (0.1)^2]);% incertidumbre inicial en bias

% Matriz de transicion F (modelo de propagacion del error)
F_lc = eye(9);
F_lc(1:3, 4:6) = eye(3) * T_s;   % pos = pos + vel*dt
F_lc(1:3, 7:9) = 0.5 * eye(3) * T_s^2; % pos = pos + 0.5*bias*dt^2
F_lc(4:6, 7:9) = eye(3) * T_s;   % vel = vel + bias*dt (error = INS - real, a_medida = a + bias)

% Ruido de proceso Q. Esto es clave para el tuning.
Q_vel = eye(3) * sigma_acc_noise^2 * T_s;
Q_bias = eye(3) * sigma_acc_drift^2 * T_s;
Q_lc = zeros(9);
Q_lc(4:6, 4:6) = Q_vel;
Q_lc(7:9, 7:9) = Q_bias;

% Matriz de observacion H. Mide el error de pos y vel.
H_lc = [eye(6), zeros(6, 3)];

% Ruido de la medida R (del GNSS)
R_lc = diag([repmat(sigma_gnss_pos^2, 1, 3), repmat(sigma_gnss_vel^2, 1, 3)]);

% Para guardar los resultados
pos_integrado_enu = zeros(num_pasos_filtro, 3);
vel_integrado_enu = zeros(num_pasos_filtro, 3);
bias_estimado_hist = zeros(num_pasos_filtro, 3);
P_diag_hist_lc = zeros(9, num_pasos_filtro);

% Bucle del filtro
for k = 1:num_pasos_filtro
    idx_real = idx_muestreo(k);
    
    pos_ins_k = pos_ins_prop(idx_real, :)';
    vel_ins_k = vel_ins_prop(idx_real, :)';
    
    % --- Prediccion ---
    x_pred = F_lc * x_lc;
    P_pred = F_lc * P_lc * F_lc' + Q_lc;
    
    % --- Actualizacion ---
    % La "medida" es la diferencia INS - GNSS
    z_medida = [pos_ins_k - pos_gnss_enu(k, :)';
                vel_ins_k - vel_gnss_enu(k, :)'];
    
    y = z_medida - H_lc * x_pred;
    
    S = H_lc * P_pred * H_lc' + R_lc;
    K = P_pred * H_lc' / S; % Ganancia de Kalman
    
    x_lc = x_pred + K * y;
    P_lc = (eye(9) - K * H_lc) * P_pred;
    
    % --- Correccion de la navegacion ---
    pos_integrado_enu(k, :) = (pos_ins_k - x_lc(1:3))';
    vel_integrado_enu(k, :) = (vel_ins_k - x_lc(4:6))';
    
    % Corregimos la trayectoria del INS para el futuro (Feedback)
    % La correccion de velocidad tambien afecta a la posicion futura: delta_v*(t - t_k)
    if idx_real < size(pos_ins_prop, 1)
        num_puntos_restantes = size(pos_ins_prop, 1) - idx_real + 1;
        dt_restante = t_real(idx_real:end) - t_real(idx_real);
        pos_ins_prop(idx_real:end, :) = pos_ins_prop(idx_real:end, :) - repmat(x_lc(1:3)', num_puntos_restantes, 1) ...
                                        - dt_restante * x_lc(4:6)';
        vel_ins_prop(idx_real:end, :) = vel_ins_prop(idx_real:end, :) - repmat(x_lc(4:6)', num_puntos_restantes, 1);
    end

    % Reseteamos los errores de pos y vel porque ya los hemos corregido
    x_lc(1:6) = 0;
    
    % Guardamos el bias estimado para verlo luego
    bias_estimado_hist(k, :) = x_lc(7:9)';
    P_diag_hist_lc(:, k) = diag(P_lc);
end

fprintf('Terminado.\n');

%% 5. Graficos para comparar
fprintf('Generando graficos...\n');

pos_real_muestreada = pos_enu_real(idx_muestreo, :);
vel_real_muestreada = vel_enu_real(idx_muestreo, :);
bias_real_muestreado = bias_acc_real(idx_muestreo, :);

% Trayectoria 2D
figure;
plot(pos_real_muestreada(:,1), pos_real_muestreada(:,2), 'k', 'LineWidth', 3);
hold on;
plot(pos_gnss_enu(:,1), pos_gnss_enu(:,2), 'm.', 'MarkerSize', 10);
plot(pos_ins_only_enu(idx_muestreo,1), pos_ins_only_enu(idx_muestreo,2), 'b--', 'LineWidth', 1.5);
plot(pos_integrado_enu(:,1), pos_integrado_enu(:,2), 'r-', 'LineWidth', 2);
grid on; axis equal;
xlabel('Posición Este [m]'); ylabel('Posición Norte [m]');
title('Comparativa de Trayectorias (Este vs. Norte)');
legend('Verdad Terreno', 'Solo GNSS (ruidosa)', 'Solo INS (con deriva)', 'Solución Integrada', 'Location', 'best');
sgtitle('Resultado de la Integración Loosely Coupled');

% Errores de posicion
error_pos_ins = pos_ins_only_enu(idx_muestreo,:) - pos_real_muestreada;
error_pos_gnss = pos_gnss_enu - pos_real_muestreada;
error_pos_integrado = pos_integrado_enu - pos_real_muestreada;

figure;
subplot(3,1,1);
plot(t_filtro, error_pos_ins(:,1), 'b--', t_filtro, error_pos_gnss(:,1), 'm:', t_filtro, error_pos_integrado(:,1), 'r-','LineWidth',1.5);
title('Error de Posición en EJE ESTE'); ylabel('[m]'); grid on;
legend('Solo INS', 'Solo GNSS', 'Integrado');
subplot(3,1,2);
plot(t_filtro, error_pos_ins(:,2), 'b--', t_filtro, error_pos_gnss(:,2), 'm:', t_filtro, error_pos_integrado(:,2), 'r-','LineWidth',1.5);
title('Error de Posición en EJE NORTE'); ylabel('[m]'); grid on;
subplot(3,1,3);
plot(t_filtro, error_pos_ins(:,3), 'b--', t_filtro, error_pos_gnss(:,3), 'm:', t_filtro, error_pos_integrado(:,3), 'r-','LineWidth',1.5);
title('Error de Posición en EJE ARRIBA'); ylabel('[m]'); xlabel('Tiempo [s]'); grid on;
sgtitle('Comparativa de Errores de Posición');

% Errores de velocidad
error_vel_ins = vel_ins_only_enu(idx_muestreo,:) - vel_real_muestreada;
error_vel_gnss = vel_gnss_enu - vel_real_muestreada;
error_vel_integrado = vel_integrado_enu - vel_real_muestreada;

figure;
subplot(3,1,1);
plot(t_filtro, error_vel_ins(:,1), 'b--', t_filtro, error_vel_gnss(:,1), 'm:', t_filtro, error_vel_integrado(:,1), 'r-','LineWidth',1.5);
title('Error de Velocidad en EJE ESTE'); ylabel('[m/s]'); grid on;
legend('Solo INS', 'Solo GNSS', 'Integrado');
subplot(3,1,2);
plot(t_filtro, error_vel_ins(:,2), 'b--', t_filtro, error_vel_gnss(:,2), 'm:', t_filtro, error_vel_integrado(:,2), 'r-','LineWidth',1.5);
title('Error de Velocidad en EJE NORTE'); ylabel('[m/s]'); grid on;
subplot(3,1,3);
plot(t_filtro, error_vel_ins(:,3), 'b--', t_filtro, error_vel_gnss(:,3), 'm:', t_filtro, error_vel_integrado(:,3), 'r-','LineWidth',1.5);
title('Error de Velocidad en EJE ARRIBA'); ylabel('[m/s]'); xlabel('Tiempo [s]'); grid on;
sgtitle('Comparativa de Errores de Velocidad');

% Estimacion del bias
figure;
sigma_bias = sqrt(P_diag_hist_lc(7:9,:));
for i = 1:3
    subplot(3,1,i)
    plot(t_filtro, bias_real_muestreado(:,i), 'k', 'LineWidth', 2);
    hold on;
    plot(t_filtro, bias_estimado_hist(:,i), 'r--','LineWidth', 1.5);
    plot(t_filtro, bias_estimado_hist(:,i) + sigma_bias(i,:)', 'b:');
    plot(t_filtro, bias_estimado_hist(:,i) - sigma_bias(i,:)', 'b:');
    grid on;
    if i==1, title('Bias Acelerómetro EJE ESTE'); end
    if i==2, title('Bias Acelerómetro EJE NORTE'); end
    if i==3, title('Bias Acelerómetro EJE ARRIBA'); xlabel('Tiempo [s]'); end
    ylabel('[m/s^2]');
    legend('Bias Real', 'Bias Estimado', '\pm1\sigma');
end
sgtitle('Seguimiento del Bias del Acelerómetro');