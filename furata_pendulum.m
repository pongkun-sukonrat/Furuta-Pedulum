% 3D Furuta Pendulum Simulation (A4 Footprint)
% Angle Convention: 90 deg = UP | 0 deg = Horizontal | -90 deg = DOWN
clear; clc; close all;

% --- 1. System Parameters ---
p.m_p = 0.1;                   % Pendulum mass (kg)
% Scaled down to fit within an A4 area (210mm x 297mm)
p.L_r = 0.10;                   % Rotary arm length (m) - 20cm diameter fits A4
p.L_p = 0.10;                   % Pendulum length (m) 
p.l_p = p.L_p / 2;              % Pendulum COM (m)
p.J_r = 0.002;                  % Arm inertia
p.J_p = (1/3) * p.m_p * (p.L_p)^2; % Pendulum inertia
p.g   = 9.81;                   % Gravity

p.J_total_arm  = p.J_r + p.m_p * p.L_r^2;
p.J_total_pend = p.J_p + p.m_p * p.l_p^2;
p.M0           = p.m_p * p.L_r * p.l_p;

% Friction & Motor Parameters
p.b_theta = 0.005; p.b_alpha = 0.001;
p.R_m     = 2.6;   p.K_t     = 0.00767;
p.K_m     = 0.00767;
p.N       = 15.0;
p.V_max   = 24.0;

p.b_eq = p.b_theta + (p.N^2 * p.K_t * p.K_m) / p.R_m;
p.gamma = (p.N * p.K_t) / p.R_m;

% --- 2. Linearized Model & LQR Gain ---
M_lin = [p.J_total_arm, p.M0; p.M0, p.J_total_pend];
invM0 = inv(M_lin);

A = zeros(4,4);
A(1,3) = 1.0; A(2,4) = 1.0;
A(3:4, 2)   = invM0 * [0; p.m_p * p.g * p.l_p];
A(3:4, 3:4) = invM0 * [-p.b_eq, 0; 0, -p.b_alpha];

B = zeros(4,1);
B(3:4, 1) = invM0 * [p.gamma; 0];

% FIXED LQR: Balanced to catch the pendulum successfully but settle gently
Q = diag([100, 600, 20, 20]); 
R = 2.0; 

K_lqr = lqr(A, B, Q, R);

% --- 3. Run Simulation ---
x0 = [0; -pi/2; 0; 0.05];
tspan = [0 10];

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
[t, x] = ode45(@(t, x) furuta_dynamics(t, x, p, K_lqr), tspan, x0, options);

% --- 4. Plot 2D State Histories ---
figure('Name', 'System Response', 'Position', [50, 40, 550, 600]);

theta_plot = rad2deg(mod(x(:,1) + pi, 2*pi) - pi);
alpha_plot = rad2deg(mod(x(:,2) + pi, 2*pi) - pi);

subplot(3,1,1); plot(t, theta_plot, 'b', 'LineWidth', 1.5);
ylabel('Arm (deg)'); grid on; title('Rotary Arm Position (\theta)');
yline(0, 'k--', 'Target (0°)'); ylim([-180 180]); yticks([-180 -90 0 90 180]);

subplot(3,1,2); plot(t, alpha_plot, 'r', 'LineWidth', 1.5);
ylabel('Pend (deg)'); grid on; title('Pendulum Position (\alpha)');
yline(90, 'k--', 'Target (90°)'); ylim([-180 180]); yticks([-180 -90 0 90 180]);

Vm = zeros(length(t), 1);
for i = 1:length(t)
    theta_user = x(i,1); alpha_user = x(i,2); dtheta = x(i,3); dalpha = x(i,4);
    alpha_std = alpha_user - pi/2;
    alpha_wrap = mod(alpha_std + pi, 2*pi) - pi;
    
    V_lqr = -K_lqr * [theta_user; alpha_wrap; dtheta; dalpha];
    E = 0.5 * p.J_total_pend * dalpha^2 + p.m_p * p.g * p.l_p * (cos(alpha_std) - 1);
    
    % FIXED SWING-UP: Massively increased energy gain and reduced arm centering constraint
    V_swing = 4500 * E * tanh(15 * dalpha * cos(alpha_std)) - 1.5 * theta_user - 1.0 * dtheta;
    
    w = exp(-(alpha_wrap / 0.5)^2);
    V_m_raw = w * V_lqr + (1 - w) * V_swing;
    Vm(i) = p.V_max * tanh(V_m_raw / p.V_max);
end

subplot(3,1,3); plot(t, Vm, 'k', 'LineWidth', 1.5);
ylabel('Voltage (V)'); xlabel('Time (s)'); grid on;

% --- 5. True 3D Animation Plot ---
figure('Name', '3D Furuta Pendulum Animation', 'Position', [620, 100, 650, 600]);

circle_angles = linspace(0, 2*pi, 100);
circle_x = p.L_r * cos(circle_angles);
circle_y = p.L_r * sin(circle_angles);

for i = 1:5:length(t)
    theta_i = x(i, 1);
    alpha_i = x(i, 2);
    
    arm_x = p.L_r * cos(theta_i);
    arm_y = p.L_r * sin(theta_i);
    arm_z = 0;
    
    pend_x = arm_x + p.L_p * cos(alpha_i) * sin(theta_i);
    pend_y = arm_y - p.L_p * cos(alpha_i) * cos(theta_i);
    pend_z = arm_z + p.L_p * sin(alpha_i);
    
    plot3([0 0], [0 0], [-0.15 0], 'k-', 'LineWidth', 6); hold on;
    plot3(circle_x, circle_y, zeros(size(circle_x)), 'k--', 'LineWidth', 0.8);
    
    plot3([0, arm_x], [0, arm_y], [0, arm_z], 'b-', 'LineWidth', 5);
    plot3([arm_x, pend_x], [arm_y, pend_y], [arm_z, pend_z], 'r-', 'LineWidth', 3.5);
    
    plot3(0, 0, -0.15, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.2 0.2 0.2]);
    plot3(0, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot3(arm_x, arm_y, arm_z, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot3(pend_x, pend_y, pend_z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Draw an A4 bounding box approximation (Green outline on the floor)
    plot3([-0.105 0.105 0.105 -0.105 -0.105], [-0.1485 -0.1485 0.1485 0.1485 -0.1485], [-0.15 -0.15 -0.15 -0.15 -0.15], 'g--', 'LineWidth', 1.5);
    hold off;
    
    axis equal;
    axis([-0.25 0.25 -0.25 0.25 -0.2 0.2]);
    grid on; box on;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    view(45, 25);
    
    theta_disp = mod(theta_i + pi, 2*pi) - pi;
    alpha_disp = mod(alpha_i + pi, 2*pi) - pi;
    
    title(sprintf('Time: %.2f s | Arm: %.1f° | Pendulum: %.1f°\nGreen Box = A4 Area', ...
        t(i), rad2deg(theta_disp), rad2deg(alpha_disp)));
    drawnow;
end

% --- Dynamics & Control Function ---
function dxdt = furuta_dynamics(~, x, p, K_lqr)
    theta_user = x(1); alpha_user = x(2); dtheta = x(3); dalpha = x(4);
    
    alpha_std = alpha_user - pi/2;
    alpha_wrap = mod(alpha_std + pi, 2*pi) - pi;
    
    V_lqr = -K_lqr * [theta_user; alpha_wrap; dtheta; dalpha];
    
    E = 0.5 * p.J_total_pend * dalpha^2 + p.m_p * p.g * p.l_p * (cos(alpha_std) - 1);
    
    % FIXED SWING-UP: Scaled up gain and loosened arm centering
    V_swing = 4500 * E * tanh(15 * dalpha * cos(alpha_std)) - 1.5 * theta_user - 1.0 * dtheta;
    
    w = exp(-(alpha_wrap / 0.5)^2);
    V_m_raw = w * V_lqr + (1 - w) * V_swing;
    
    V_m = p.V_max * tanh(V_m_raw / p.V_max);
    
    tau_arm = p.gamma * V_m - p.b_eq * dtheta;
    
    M11 = p.J_total_arm + p.m_p * p.l_p^2 * (sin(alpha_std))^2;
    M12 = p.m_p * p.L_r * p.l_p * cos(alpha_std);
    M = [M11, M12; M12, p.J_total_pend];
    
    C1 = p.m_p * p.l_p^2 * sin(2*alpha_std) * dtheta * dalpha - p.m_p * p.L_r * p.l_p * sin(alpha_std) * dalpha^2;
    C2 = -0.5 * p.m_p * p.l_p^2 * sin(2*alpha_std) * dtheta^2 - p.m_p * p.g * p.l_p * sin(alpha_std) + p.b_alpha * dalpha;
    
    q_dd = M \ [tau_arm - C1; -C2];
    dxdt = [dtheta; dalpha; q_dd(1); q_dd(2)];
end








