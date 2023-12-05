function [lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T)

% phi_i: unmet demand penalty
phi_l = 20;
phi_h = 30;
phi = phi_l + (phi_h - phi_l).* rand(I,1); 

% demand at each time slot
lambda_l = 3;
lambda_h = 10;
lambda_fore = round(lambda_l + (lambda_h - lambda_l).*rand(I,T));

% Number of data center: C_max(j)
C_max = zeros(J,1);
C_max_set = [50 100 200];
for j = 1:J
    % Select the capacity based on the population
    index = randi(3,1);
    C_max(j) = C_max_set(index);
end 


% average server utilization control 
gamma_max = 0.9;
rho = 0.8;

% Time interval
Delta_T = 1;

% idle power: ~ 0.45 kwh
P_idle = 0.45 + (0.55 - 0.45) * rand(J,1);
% P_idle = P_idle * 10;
% computation: 1.2 kwh
P_peak = 1.2 + (1.5 - 1.2) * rand(J,1); 
% P_peak = P_peak * 10;
% PUE: Energy efficiency
E_usage = 1.8 + (1.9 - 1.8) * rand(J,1);

% CSP maximum grid power
P_G_h = 1000;
P_G_l = 800;
P_G_max = P_G_l + (P_G_h - P_G_l) * rand(J,T);

% renewable energy: 
P_R_h = 60;
P_R_l = 50;
P_R = P_R_l + (P_R_h - P_R_l) * rand(J,T);

% Battery charaged and discharg rate
P_C_h = 80;
P_C_l = 70;
P_C_max = P_C_l + (P_C_h - P_C_l) * rand(J,1);

% 
P_D_h = 80;
P_D_l = 70;
P_D_max = P_D_l + (P_D_h - P_D_l) * rand(J,1);


% electricity price
e_l = 0.1;
e_h = 0.35;
e = e_l + (e_h - e_l) * rand(J,T);

% sell-back price to the grid
a = 0.8 * e;

eta_c = 0.8;
eta_d = 0.8;

% capacity of battery at DC j
E_max = 100;
E_min = 30;
E_0 = 20 * ones(J,1);

% Emission factor 

theta_grid_h = 0.8;
theta_grid_l = 0.1;
theta_grid = theta_grid_l + (theta_grid_h - theta_grid_l) * rand(J,1);

% Carbon tax: To address environmental concerns: delta_j is defined as
% carbon tax for each location j
delta_h = 0.7;
delta_l = 0.6;
delta = delta_l + (delta_h - delta_l) * rand(J,1);


% alpha: average computing resource per request
alpha = 0.1;

