% Flow over a bump in a channel 
% Written by Kevin Morales

clear all
close all
clc
tic

%_____________Domain______________
% Load mesh
load("mesh_bumpchannel2.mat");
% load("mesh_bumpchannel353_161.mat");
% Y coordinates in Y
% X coordinates in X
[n_y, n_x] = size(X); % number of cells in X and Y
n_x = n_x - 1;
n_y = n_y - 1;
dx = X(1, 2) - X(1, 1); % dx  
L_y = max(Y(:, 1));

%_____________Constants______________
Vel_initial = 1e-4;  % 29e-4; % Velocity at the inlet 
rho = 1.2; % Density (Kg/m3)
nu = 0.0000174; % Viscosity (Kg/(m*s))
Reynolds_num = Vel_initial * rho * L_y / nu; % Dimensionless 
mass_flow = L_y * Vel_initial;
disp(Reynolds_num)

%_____________Velocity and Pressure Fields ______________
U_vel = zeros(n_y + 2, n_x + 1); % Velocity field in X (m/s)
d_e = zeros(size(U_vel)); % Pressure correction coefficients for the velocity field in X
% Defining initial and boundary conditions for U
U_vel(:, 1) = Vel_initial; % Inlet velocity 
U_vel(1, :) = -U_vel(2, :) + 2 * Vel_initial; % Free stream velocity
U_vel(n_y + 2, :) = -U_vel(n_y + 1, :); % No-slip condition at the south wall

V_vel = zeros(n_y + 1, n_x + 2); % Velocity field in Y (m/s)
d_n = zeros(size(V_vel)); % Pressure correction coefficients for the velocity field in Y
% Defining initial and boundary conditions for V
V_vel(1, :) = 0; % No-slip condition at the north wall
V_vel(n_y + 1, :) = 0; % No-slip condition at the south wall

P_press = ones(n_y + 2, n_x + 2); % Pressure field
P_corr = zeros(size(P_press)); % Pressure correction field
% Defining initial and boundary conditions for P
P_press(n_y + 2, :) = P_press(n_y + 1, :); % NO mass flow at the bottom face
% P_press(1, :) = P_press(2, :); % No mass flow at the upper face
P_press(:, 1) = P_press(:, 2) + 0.25 * (U_vel(:, 1).^2); % Bernoulli for flow at the inlet

%_____________Additional Coefficients ______________
% Under-relaxation Factors
alpha_p = 0.7;
alpha_vel = 0.5;

% Error measurement through continuity equation
b_p = zeros(size(P_press));

% Initializing error
err = 1;
% Required error
err_req = 1e-7;
error_vec = [];
% Iterations
iterations = 0;

%_____________ Main Iteration  ______________
while err > err_req
    if err > 1e12
        fprintf("Unstable solution, iteration has been stopped")
        break
    end
    % Momentum equations in X, Y, and pressure correction
    for i = 1:n_y
        for j = 1:n_x - 1
            % Here go the functions
            x = [X(n_y + 1 - i, j) X(n_y + 1 - i, j + 1) X(n_y + 1 - i, j + 2) X(n_y + 2 - i, j + 2) X(n_y + 2 - i, j + 1) X(n_y + 2 - i, j)];
            y = [Y(n_y + 1 - i, j) Y(n_y + 1 - i, j + 1) Y(n_y + 1 - i, j + 2) Y(n_y + 2 - i, j + 2) Y(n_y + 2 - i, j + 1) Y(n_y + 2 - i, j)];
            m = n_y + 2 - i;
            % m = i + 1;
            n = j + 1;       
            u_N = U_vel(m - 1, n);
            u_E = U_vel(m, n + 1);
            u_S = U_vel(m + 1, n);
            u_W = U_vel(m, n - 1);
            u_P = U_vel(m, n); 
            u_N_w = U_vel(m - 1, n - 1);
            u_N_e = U_vel(m - 1, n + 1); 
            u_S_w = U_vel(m + 1, n - 1);
            u_S_e = U_vel(m + 1, n + 1);
            v_n_w = V_vel(m - 1, n);
            v_n_e = V_vel(m - 1, n + 1);
            v_s_w = V_vel(m, n);
            v_s_e = V_vel(m, n + 1);
            p_P = P_press(m, n);
            p_E = P_press(m, n + 1);
            p_P_n = P_press(m - 1, n);
            p_E_n = P_press(m - 1, n + 1);
            p_P_s = P_press(m + 1, n);
            p_E_s = P_press(m + 1, n + 1);
            
            % Solution of the momentum equation in X
            [U_vel(m, n), d_e(m, n)] = x_Momentum_U(x, y, u_N, u_E, u_S, u_W, u_P, u_N_w, u_N_e, u_S_w, u_S_e, v_n_w, v_n_e, v_s_w, v_s_e, rho, nu, p_P, p_E, p_P_n, p_E_n, p_P_s, p_E_s, dx);
            
        end
    end
    % Reapply boundary conditions in X
    U_vel(1, :) = -U_vel(2, :) + 2 * Vel_initial; % Free stream velocity
    % U_vel(n_y + 2, :) = -U_vel(n_y + 1, :); % No-slip condition at the south wall
    U_vel(:, n_x + 1) = U_vel(:, n_x); % Zero velocity gradient at the outlet

    for i = 1:n_y - 1
        for j = 1:n_x
            x2 = [X(n_y - i, j) X(n_y - i, j + 1) X(n_y + 1 - i, j + 1) X(n_y + 2 - i, j + 1) X(n_y + 2 - i, j) X(n_y + 1 - i, j)];
            y2 = [Y(n_y - i, j) Y(n_y - i, j + 1) Y(n_y + 1 - i, j + 1) Y(n_y + 2 - i, j + 1) Y(n_y + 2 - i, j) Y(n_y + 1 - i, j)];
            
            m = n_y + 1 - i;
            % m = i + 1;
            n = j + 1;
            v_N = V_vel(m - 1, n);
            v_E = V_vel(m, n + 1);
            v_S = V_vel(m + 1, n);
            v_W = V_vel(m, n - 1);
            v_P = V_vel(m, n); 
            v_N_w = V_vel(m - 1, n - 1);
            v_N_e = V_vel(m - 1, n + 1);
            v_S_w = V_vel(m + 1, n - 1);
            v_S_e = V_vel(m + 1, n - 1);
            u_n_w = U_vel(m, n - 1);
            u_n_e = U_vel(m, n);
            u_s_w = U_vel(m + 1, n - 1);
            u_s_e = U_vel(m + 1, n);
            p_P = P_press(m, n);
            p_S = P_press(m + 1, n);
            % Solution of the momentum equation in Y
            [V_vel(m, n), d_n(m, n)] = y_Momentum_V(x2, y2, v_N, v_E, v_S, v_W, v_P, v_N_w, v_N_e, v_S_w, v_S_e, u_n_w, u_n_e, u_s_w, u_s_e, rho, nu, p_P, p_S, dx);
            
        end 
    end
    % Boundary conditions for Y
    V_vel(:, 1) = -V_vel(:, 2); % No Y velocity component at the inlet
    % V_vel(:, n_x + 2) = V_vel(:, n_x + 1); % Zero velocity gradient at the outlet
    % V_vel(1, :) = 0; % No mass flow at the upper face
    V_vel(n_y + 1, :) = 0; % No-slip condition at the south wall

    % Pressure Correction 
    P_corr = zeros(size(P_press));
    for i = 1:n_y
        for j = 1:n_x
            x3 = [X(n_y + 1 - i, j), X(n_y + 1 - i, j + 1), X(n_y + 2 - i, j + 1), X(n_y + 2 - i, j)];
            y3 = [Y(n_y + 1 - i, j), Y(n_y + 1 - i, j + 1), Y(n_y + 2 - i, j + 1), Y(n_y + 2 - i, j)];

            m = n_y + 2 - i;
            n = j + 1;
            U_e = U_vel(m, n);
            U_w = U_vel(m, n - 1);
            V_n = V_vel(m - 1, n);
            V_s = V_vel(m, n);
            u_N_w = U_vel(m - 1, n - 1);
            u_N_e = U_vel(m - 1, n);
            u_S_w = U_vel(m + 1, n - 1);
            u_S_e = U_vel(m + 1, n);
            P_n = P_corr(m - 1, n);
            P_e = P_corr(m, n + 1);
            P_s = P_corr(m + 1, n);
            P_w = P_corr(m, n - 1);
            a_e = -d_e(m, n); 
            a_w = -d_e(m, n - 1);
            a_s = -d_n(m, n);
            a_n = -d_n(m - 1, n);
            [P_corr(m, n), b_p(m, n)] = Pressure_corr(x3, y3, a_e, a_w, a_s, a_n, U_e, U_w, V_n, V_s, u_N_w, u_N_e, u_S_w, u_S_e, P_n, P_e, P_s, P_w);
        end 
    end
    P_press = P_press + alpha_p * P_corr;
    % Reapply boundary conditions to the pressure
    P_press(:, 1) = P_press(:, 2) + 0.25 * (U_vel(:, 1).^2); % Pressure difference required for velocity condition at the flow inlet
    P_press(n_y + 2, :) = P_press(n_y + 1, :); % NO mass flow at the bottom face
    % P_press(1, :) = P_press(2, :); % No mass flow at the upper face

    % Correction of velocity in X
    for i = 1:n_y
        for j = 1:n_x - 1
            m = n_y + 2 - i;
            n = j + 1;
            U_vel(m, n) = U_vel(m, n) + alpha_vel * d_e(m, n) * (P_corr(m, n + 1) - P_corr(m, n));
        end
    end
    % Reapply boundary conditions in U
    U_vel(1, :) = -U_vel(2, :) + 2 * Vel_initial; % Free stream velocity
    U_vel(n_y + 2, :) = -U_vel(n_y + 1, :); % No-slip condition at the south wall
    U_vel(:, n_x + 1) = U_vel(:, n_x); % Zero velocity gradient at the outlet

    % Correction of velocity in Y
    for i = 1:n_y - 1
        for j = 1:n_x
            m = n_y + 1 - i;
            % m = i + 1;
            n = j + 1;
            V_vel(m, n) = V_vel(m, n) + alpha_vel * d_n(m, n) * (P_corr(m, n) - P_corr(m + 1, n));
        end
    end
    % Reapply boundary conditions in V
    % Reapply boundary conditions in V
    V_vel(:, 1) = -V_vel(:, 2); % No Y velocity component at the inlet
    V_vel(1, :) = 0; % No-slip condition at the north wall
    V_vel(n_y + 1, :) = 0; % No-slip condition at the south wall

    err = 0;
    % Calculate the error as the residual of the continuity equation
    for i = 1:n_y + 2
        for j = 1:n_x + 2
            err = err + abs(b_p(i, j));
        end
    end
    err;
    iterations = iterations + 1
    error_vec(iterations) = err;
    
end 
beep();
toc

save('bumpinchannel_fields2.mat', 'U_vel', 'V_vel', 'P_press', 'b_p');

iterations_vec = 1:iterations;

figure(26);
plot(iterations_vec, error_vec);
title("Error");
xlabel("Number of iterations");
ylabel("Error");
