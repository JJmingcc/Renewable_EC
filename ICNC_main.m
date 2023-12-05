% Simulation setup 

I = 10; 
J = 8;
T = 12;



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
delta_l = 0.6;
delta = delta_l + (delta_h - delta_l) * rand(J,1);


% alpha: average computing resource per request
alpha = 0.1;
delta_h = 0.7;


%% Initialization

I = 10; 
J = 8;
T = 12;
 [lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);



%% Cloud service provider optimization model

cvx_begin
    cvx_solver Gurobi_2
    cvx_precision low
    % computing installation
    variable c(J,T) integer
    % workload allocation
    variable x(I,J,T)
    variable q(I,T)
    variable workload(J,T)
    % utilization
    variable gamma_var(J,T)
    % power input
    variable P_G(J,T)
    variable P_C(J,T)
    variable P_D(J,T)
    variable P_S(J,T)
    variable P_U(J,T)
    variable EM(J,T)
    variable E(J,T)
    % master obj funtion
    variable C_u
    variable C_e
    variable C_c
    variable obj
    minimize obj
    subject to 
          obj == C_u + C_c + C_e;
            % Demand supply balance: Eq(2)
            for t = 1:T
                alpha * sum(x(:,:,t),2) + q(:,t) == lambda_fore(:,t); 
            end
            % Obj1: C^u: Eq(3)
            C_u == sum(sum((phi * ones(1,T)).*q));
            % Obj2: Battery recharge and discharge Eq(9)
            C_e == sum(sum(e .* P_G - a.* P_S));      
            gamma_var == reshape(sum(x,1),[J,T]);
            
            % Obj3: Carbon tax
            for j  = 1:J
                for t = 1:T
                    % Eq(1)
                    c(j,t) <= C_max(j);
                    % Eq(5)
                    gamma_var(j,t) <= rho * c(j,t) * gamma_max;
                    % Eq(6)
                    P_U(j,t) == c(j,t) * (P_idle(j) + (E_usage(j) - 1)* P_peak(j)) + (P_idle(j) - P_peak(j)) * gamma_var(j,t)/ rho; % c cancel out
                    % Eq(14)
                    EM(j,t) == theta_grid(j) * P_G(j,t);
                end
            end   
            C_c == sum(sum((delta * ones(1,T)) .* EM));
            % Power blance equation
            for t = 1:T
                for j = 1:J
                    P_G(j,t) + P_R(j,t) + P_D(j,t) == P_U(j,t) + P_C(j,t) + P_S(j,t);
                    P_C(j,t) <= P_C_max(j);
                    P_D(j,t) <= P_D_max(j);
                    P_G(j,t) <= P_G_max(j,t);
                end
            end

            % Eq(12) & (13)
            E <= E_max;
            E >= E_min;
            for t = 1:T-1
                for j  = 1:J
                    if t == 1
                       E(j,1) == E_0(j,1) + Delta_T * (eta_c * P_C(j,t) - P_D(j,t)/eta_d);
                    else 
                       E(j,t+1) == E(j,t) + Delta_T * (eta_c * P_C(j,t) - P_D(j,t)/eta_d);
                    end
                end
            end
            c >= 0;
            x >= 0;
            q >= 0;
            P_G >= 0;
            P_C >= 0;
            P_D >= 0;
            P_S >= 0;
cvx_end

%% Experiment

[cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);

%% Sensitivity analysis on renewable energy and price 

I = 10; 
J = 8;
T = 12;
Psi_set = 0.2:0.2:1.6;
xi_e_set = 0.4:0.2:1.4;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);

Psi_e_cost_set = zeros(length(Psi_set),length(xi_e_set));
P_R_fix = P_R;
e_fix = e;
for i = 1:length(Psi_set)
     Psi = Psi_set(i);
     P_R = P_R_fix * Psi;
     for j = 1: length(xi_e_set)
        xi_e = xi_e_set(j);
        e = xi_e * e_fix;
        [Psi_e_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
        Psi_e_cost_set(i,j) =  Psi_e_cost;
     end
end


%%

Psi_set = 0.2:0.2:1.6;
figure;
plot(Psi_set,Psi_e_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(Psi_set,Psi_e_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(Psi_set,Psi_e_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(Psi_set,Psi_e_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(Psi_set,Psi_e_cost_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(Psi_set,Psi_e_cost_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks([0.2:0.2:1.6]);
xlabel('\Psi','FontSize',18);
legend('\xi_e = 0.4','\xi_e = 0.6','\xi_e = 0.8','\xi_e = 1.0','\xi_e = 1.2','\xi_e = 1.4','FontSize',18)
ylabel('Total cost','FontSize',18);


%% Sensitivity analysis on E_max and price 

I = 10; 
J = 8;
T = 12;
E_max_ratio_set = 0.6:0.2:1.8;
xi_e_set = 0.4:0.2:1.4;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
E_e_cost_set = zeros(length(E_max_ratio_set),length(xi_e_set));

E_max_fix = E_max;
e_fix = e;

for i = 1:length(E_max_ratio_set)
     E_max_ratio = E_max_ratio_set(i);
     E_max = E_max_fix * E_max_ratio;
     for j = 1: length(xi_e_set)
        xi_e = xi_e_set(j);
        e = xi_e * e_fix;
        [E_e_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
        E_e_cost_set(i,j) =  E_e_cost;
     end
end

%%
E_max_ratio_set = 0.6:0.2:1.8;
figure;
plot(E_max_ratio_set,E_e_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,E_e_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(E_max_ratio_set,E_e_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,E_e_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(E_max_ratio_set,E_e_cost_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(E_max_ratio_set,E_e_cost_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.6 1.8]);
ylim([4800 6000]);

xticks([0.6:0.2:1.8]);
xlabel('\xi_{E^{max}}','FontSize',18);
legend('\xi_e = 0.4','\xi_e = 0.6','\xi_e = 0.8','\xi_e = 1.0','\xi_e = 1.2','\xi_e = 1.4','FontSize',16)
ylabel('Total cost','FontSize',18);


%% Sensitivity analysis on E_max and price 

I = 10; 
J = 8;
T = 12;
Psi_set = 0.2:0.2:1.6;
zeta_set = 0:0.2:1.0;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);

zeta_Psi_rev_set = zeros(length(Psi_set),length(zeta_set));
zeta_Psi_cost_set = zeros(length(Psi_set),length(zeta_set));

for i = 1:length(Psi_set)
     Psi = Psi_set(i);
     P_R = P_R_fix * Psi;
     for j = 1: length(zeta_set)
        zeta = zeta_set(j);
        a = e * zeta;
        [zeta_Psi_rev,zeta_Psi_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
        zeta_Psi_rev_set(i,j) =  zeta_Psi_rev;
        zeta_Psi_cost_set(i,j) = zeta_Psi_cost;
     end
end

%% 

figure;
plot(Psi_set,zeta_Psi_rev_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(Psi_set,zeta_Psi_rev_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(Psi_set,zeta_Psi_rev_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(Psi_set,zeta_Psi_rev_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(Psi_set,zeta_Psi_rev_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(Psi_set,zeta_Psi_rev_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks([0.2:0.2:1.6]);
xlabel('\Psi','FontSize',18);
legend('\zeta = 0.0','\zeta = 0.2','\zeta = 0.4','\zeta = 0.6','\zeta = 0.8','\zeta = 1.0','FontSize',18)
ylabel('Total cost','FontSize',18);


%%

figure;
plot(Psi_set,zeta_Psi_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(Psi_set,zeta_Psi_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(Psi_set,zeta_Psi_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(Psi_set,zeta_Psi_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(Psi_set,zeta_Psi_cost_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(Psi_set,zeta_Psi_cost_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);
xlim([0.2 1.6]);
xticks([0.2:0.2:1.6]);
xlabel('\Psi','FontSize',18);
legend('\zeta = 0.0','\zeta = 0.2','\zeta = 0.4','\zeta = 0.6','\zeta = 0.8','\zeta = 1.0','FontSize',18)
ylabel('Total cost','FontSize',18);

%% Sensitivity analysis on E_max and price 

I = 10; 
J = 8;
T = 12;
zeta_set = 0:0.2:1.0;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
E_max_ratio_set = 0.6:0.2:1.8;
E_max_fix = E_max;
zeta_E_rev_set = zeros(length(E_max_ratio_set),length(zeta_set));
zeta_E_cost_set = zeros(length(E_max_ratio_set),length(zeta_set));

for i = 1:length(E_max_ratio_set)
     E_max_ratio = E_max_ratio_set(i);
     E_max = E_max_fix * E_max_ratio;
     for j = 1: length(zeta_set)
        zeta = zeta_set(j);
        a = e * zeta;
        [zeta_E_rev,zeta_E_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
        zeta_E_rev_set(i,j) =  zeta_E_rev;
        zeta_E_cost_set(i,j) = zeta_E_cost;
     end
end

%%

E_max_ratio_set = 0.6:0.2:1.8;
figure;
plot(E_max_ratio_set,zeta_E_rev_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,zeta_E_rev_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(E_max_ratio_set,zeta_E_rev_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,zeta_E_rev_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(E_max_ratio_set,zeta_E_rev_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(E_max_ratio_set,zeta_E_rev_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.6 1.8]);
ylim([0 350]);
xticks([0.6:0.2:1.8]);
xlabel('\xi_{E^{max}}','FontSize',18);
legend('\zeta = 0.0','\zeta = 0.2','\zeta = 0.4','\zeta = 0.6','\zeta = 0.8','\zeta = 1.0','FontSize',18)
ylabel('Total cost','FontSize',18);

%%
E_max_ratio_set = 0.6:0.2:1.8;
figure;
plot(E_max_ratio_set,zeta_E_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,zeta_E_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(E_max_ratio_set,zeta_E_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,zeta_E_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(E_max_ratio_set,zeta_E_cost_set(:,5),'-+','linewidth',2,'markersize',16);
hold on 
plot(E_max_ratio_set,zeta_E_cost_set(:,6),'-d','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.6 1.8]);
xticks([0.6:0.2:1.8]);
xlabel('\xi_{E^{max}}','FontSize',18);
legend('\zeta = 0.0','\zeta = 0.2','\zeta = 0.4','\zeta = 0.6','\zeta = 0.8','\zeta = 1.0','FontSize',18)
ylabel('Total cost','FontSize',18);



%% 


I = 10; 
J = 8;
T = 12;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
gamma_max_set = 0.2:0.2:1;
D_max_set = [80;64;56;48;32;24;12];
gamma_b_cost_set = zeros(length(D_max_set),length(gamma_max_set));

for i = 1:length(D_max_set)
    D_max = D_max_set(i);
    b = zeros(I,J);
    b(randperm(numel(b),D_max))=1; 
     for j = 1: length(gamma_max_set)
        gamma_max = gamma_max_set(j);
        [gamma_b_rev,gamma_b_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage,b);
        gamma_b_cost_set(i,j) = gamma_b_cost;
     end
end


%% 


b_max_set = 0.2:0.1:0.8;
figure;
plot(b_max_set,gamma_b_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(b_max_set,gamma_b_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(b_max_set,gamma_b_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(b_max_set,gamma_b_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(b_max_set,gamma_b_cost_set(:,5),'-+','linewidth',2,'markersize',16);
grid on
set(gca,'FontSize',18);

xlim([0.2 0.8]);
ylim([5000 18000]);
xticks([0.2:0.1:0.8]);
xlabel('\xi_{D^{max}}','FontSize',18);
legend('\gamma^{max} = 0.2','\gamma^{max} = 0.4','\gamma^{max} = 0.6','\gamma^{max} = 0.8','\gamma^{max} = 1.0','FontSize',16)
ylabel('Total cost','FontSize',18);


%% 

EN_set = [8;12;16;20;24];
AP_set = [10;14;18;22;26;30];
EN_AP_cost_set = zeros(length(EN_set),length(AP_set));
T = 12;

[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
b = ones(I,J);

for i = 1:length(AP_set)
    I = AP_set(i);
     for j = 1: length(EN_set)
        J = EN_set(j);
        b = ones(I,J);
        [lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
        a = zeros(J,T);
        [EN_AP_rev,EN_AP_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage,b);
        EN_AP_cost_set(i,j) = EN_AP_cost;
     end
end



%%

figure;
plot(AP_set,EN_AP_cost_set(:,1),'-x','linewidth',2,'markersize',16);
hold on
plot(AP_set,EN_AP_cost_set(:,2),'-s','linewidth',2,'markersize',12);
hold on
plot(AP_set,EN_AP_cost_set(:,3),'-*','linewidth',2,'markersize',16);
hold on
plot(AP_set,EN_AP_cost_set(:,4),'-v','linewidth',2,'markersize',12);
hold on 
plot(AP_set,EN_AP_cost_set(:,5),'-+','linewidth',2,'markersize',16);
grid on
set(gca,'FontSize',18);

xlim([10 30]);
xticks(10:4:30);
xlabel('M','FontSize',18);
legend('N = 8','N = 12','N = 16','N = 20','N = 24','FontSize',16)
ylabel('Total cost','FontSize',18);




%%


E_max_ratio_set = 0.6:0.2:1.8;
M0_E_max_compare_cost = zeros(1, length(E_max_ratio_set));
M1_E_max_compare_cost = zeros(1, length(E_max_ratio_set));
M2_E_max_compare_cost = zeros(1, length(E_max_ratio_set));
M3_E_max_compare_cost = zeros(1, length(E_max_ratio_set));

M0_E_max_compare_rev = zeros(1, length(E_max_ratio_set));
M1_E_max_compare_rev = zeros(1, length(E_max_ratio_set));
M2_E_max_compare_rev = zeros(1, length(E_max_ratio_set));
M3_E_max_compare_rev = zeros(1, length(E_max_ratio_set));

I = 10;
J = 8;
T = 12;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
b = ones(I,J);
E_max_fix = E_max;
for i = 1:length(E_max_ratio_set)
    E_max_ratio = E_max_ratio_set(i);
    E_max  = E_max_fix * E_max_ratio;
%     [M0_rev,M0_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = DET_main(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage,b);  
    [M0_rev,M0_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = Model_0(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M1_cost, P_G,c,x,q,uti_percent] = Model_1(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,gamma_max,E_usage);
    [M2_cost,P_G,P_C,P_D,c,x,q,uti_percent] = Model_2(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M3_rev, M3_cost,P_G,P_S,c,x,q,uti_percent] = Model_3(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,Delta_T,gamma_max,E_usage);

    M0_E_max_compare_cost(i) = M0_cost;
    M1_E_max_compare_cost(i) = M1_cost;
    M2_E_max_compare_cost(i) = M2_cost;
    M3_E_max_compare_cost(i) = M3_cost;

    M0_E_max_compare_rev(i) = M0_rev;
    M1_E_max_compare_rev(i) = 0;
    M2_E_max_compare_rev(i) = 0;
    M3_E_max_compare_rev(i) = M3_rev;
end
%%

E_max_ratio_set = 0.6:0.2:1.8;

figure;
plot(E_max_ratio_set,M0_E_max_compare_cost,'-x','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,M1_E_max_compare_cost,'-s','linewidth',2,'markersize',12);
hold on
plot(E_max_ratio_set,M2_E_max_compare_cost,'-*','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,M3_E_max_compare_cost,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.6 1.8]);
xticks(0.6:0.2:1.8);
xlabel('\xi_{E^{max}}','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Total cost','FontSize',18);

%%
E_max_ratio_set = 0.6:0.2:1.8;

figure;
plot(E_max_ratio_set,M0_E_max_compare_rev,'-x','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,M1_E_max_compare_rev,'-s','linewidth',2,'markersize',12);
hold on
plot(E_max_ratio_set,M2_E_max_compare_rev,'-*','linewidth',2,'markersize',16);
hold on
plot(E_max_ratio_set,M3_E_max_compare_rev,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.6 1.8]);
xticks(0.6:0.2:1.8);
xlabel('\xi_{E^{max}}','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Total cost','FontSize',18);





%%


zeta_set = 0:0.2:1.0;
M0_zeta_compare_cost = zeros(1, length(zeta_set));
M1_zeta_compare_cost = zeros(1, length(zeta_set));
M2_zeta_compare_cost = zeros(1, length(zeta_set));
M3_zeta_compare_cost = zeros(1, length(zeta_set));

M0_zeta_compare_rev = zeros(1, length(zeta_set));
M1_zeta_compare_rev = zeros(1, length(zeta_set));
M2_zeta_compare_rev = zeros(1, length(zeta_set));
M3_zeta_compare_rev = zeros(1, length(zeta_set));

I = 10;
J = 8;
T = 12;
% [lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
E_max = E_max_fix;


for i = 1:length(zeta_set)
    zeta = zeta_set(i);
    a = e * zeta;
    [M0_rev,M0_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = Model_0(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M1_cost, P_G,c,x,q,uti_percent] = Model_1(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,gamma_max,E_usage);
    [M2_cost,P_G,P_C,P_D,c,x,q,uti_percent] = Model_2(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M3_rev, M3_cost,P_G,P_S,c,x,q,uti_percent] = Model_3(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,Delta_T,gamma_max,E_usage);

    M0_zeta_compare_cost(i) = M0_cost;
    M1_zeta_compare_cost(i) = M1_cost;
    M2_zeta_compare_cost(i) = M2_cost;
    M3_zeta_compare_cost(i) = M3_cost;

    M0_zeta_compare_rev(i) = M0_rev;
    M1_zeta_compare_rev(i) = 0;
    M2_zeta_compare_rev(i) = 0;
    M3_zeta_compare_rev(i) = M3_rev;
end


%%
zeta_set = 0:0.2:1.0;

figure;
plot(zeta_set,M0_zeta_compare_cost,'-x','linewidth',2,'markersize',16);
hold on
plot(zeta_set,M1_zeta_compare_cost,'-s','linewidth',2,'markersize',12);
hold on
plot(zeta_set,M2_zeta_compare_cost,'-*','linewidth',2,'markersize',16);
hold on
plot(zeta_set,M3_zeta_compare_cost,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.0 1.0]);
xticks(0.0:0.2:1.0);
xlabel('\zeta','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Total cost','FontSize',18);

%%
zeta_set = 0:0.2:1.0;

figure;
plot(zeta_set,M0_zeta_compare_rev,'-x','linewidth',2,'markersize',16);
hold on
plot(zeta_set,M1_zeta_compare_rev,'-s','linewidth',2,'markersize',12);
hold on
plot(zeta_set,M2_zeta_compare_rev,'-*','linewidth',2,'markersize',16);
hold on
plot(zeta_set,M3_zeta_compare_rev,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.0 1.0]);
xticks(0.0:0.2:1.0);
xlabel('\zeta','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Revenue','FontSize',18);


%%
Psi_set = 0.2:0.2:1.6;
M0_Psi_compare_cost = zeros(1, length(Psi_set));
M1_Psi_compare_cost = zeros(1, length(Psi_set));
M2_Psi_compare_cost = zeros(1, length(Psi_set));
M3_Psi_compare_cost = zeros(1, length(Psi_set));

M0_Psi_compare_rev = zeros(1, length(Psi_set));
M1_Psi_compare_rev = zeros(1, length(Psi_set));
M2_Psi_compare_rev = zeros(1, length(Psi_set));
M3_Psi_compare_rev = zeros(1, length(Psi_set));

% I = 10;
% J = 8;
% T = 12;
[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
a = 0.8 * e;
Psi_fix = Psi;
for i = 1:length(Psi_set)
    Psi_ratio = Psi_set(i);
    Psi = Psi_fix * Psi_ratio;
    [M0_rev,M0_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = Model_0(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M1_cost, P_G,c,x,q,uti_percent] = Model_1(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,gamma_max,E_usage);
    [M2_cost,P_G,P_C,P_D,c,x,q,uti_percent] = Model_2(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M3_rev, M3_cost,P_G,P_S,c,x,q,uti_percent] = Model_3(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,Delta_T,gamma_max,E_usage);

    M0_Psi_compare_cost(i) = M0_cost;
    M1_Psi_compare_cost(i) = M1_cost;
    M2_Psi_compare_cost(i) = M2_cost;
    M3_Psi_compare_cost(i) = M3_cost;

    M0_Psi_compare_rev(i) = M0_rev;
    M1_Psi_compare_rev(i) = 0;
    M2_Psi_compare_rev(i) = 0;
    M3_Psi_compare_rev(i) = M3_rev;
end


%%
Psi_set = 0.2:0.2:1.6;

figure;
plot(Psi_set,M0_Psi_compare_cost,'-x','linewidth',2,'markersize',16);
hold on
plot(Psi_set,M1_Psi_compare_cost,'-s','linewidth',2,'markersize',12);
hold on
plot(Psi_set,M2_Psi_compare_cost,'-*','linewidth',2,'markersize',16);
hold on
plot(Psi_set,M3_Psi_compare_cost,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks(0.2:0.2:1.6);
xlabel('\Psi','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Total cost','FontSize',18);

%%
Psi_set = 0.2:0.2:1.6;

figure;
plot(Psi_set,M0_Psi_compare_rev,'-x','linewidth',2,'markersize',16);
hold on
plot(Psi_set,M1_Psi_compare_rev,'-s','linewidth',2,'markersize',12);
hold on
plot(Psi_set,M2_Psi_compare_rev,'-*','linewidth',2,'markersize',16);
hold on
plot(Psi_set,M3_Psi_compare_rev,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks(0.2:0.2:1.6);
xlabel('\Psi','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Revenue','FontSize',18);


%%
xi_e_set = 0.2:0.2:1.6;
M0_e_compare_cost = zeros(1, length(xi_e_set));
M1_e_compare_cost = zeros(1, length(xi_e_set));
M2_e_compare_cost = zeros(1, length(xi_e_set));
M3_e_compare_cost = zeros(1, length(xi_e_set));

M0_e_compare_rev = zeros(1, length(xi_e_set));
M1_e_compare_rev = zeros(1, length(xi_e_set));
M2_e_compare_rev = zeros(1, length(xi_e_set));
M3_e_compare_rev = zeros(1, length(xi_e_set));


[lambda_fore,C_max,Delta_T,P_idle,P_peak,E_usage,P_G_max,P_R, P_C_max,P_D_max,alpha,rho,eta_c,eta_d,e,a,theta_grid,delta,gamma_max,phi,E_0,E_max,E_min] = init(I,J,T);
e_fix  = e;

for i = 1:length(xi_e_set)
    xi_e = xi_e_set(i);
    e = e_fix * xi_e;
    [M0_rev,M0_cost,P_G,P_C,P_D,P_S,c,x,q,uti_percent] = Model_0(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M1_cost, P_G,c,x,q,uti_percent] = Model_1(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,gamma_max,E_usage);
    [M2_cost,P_G,P_C,P_D,c,x,q,uti_percent] = Model_2(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,P_idle,P_peak,P_R,P_C_max,P_G_max,P_D_max,E_max,E_min,E_0,delta,theta_grid,eta_c,eta_d,Delta_T,gamma_max,E_usage);
    [M3_rev, M3_cost,P_G,P_S,c,x,q,uti_percent] = Model_3(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,Delta_T,gamma_max,E_usage);

    M0_e_compare_cost(i) = M0_cost;
    M1_e_compare_cost(i) = M1_cost;
    M2_e_compare_cost(i) = M2_cost;
    M3_e_compare_cost(i) = M3_cost;

    M0_e_compare_rev(i) = M0_rev;
    M1_e_compare_rev(i) = 0;
    M2_e_compare_rev(i) = 0;
    M3_e_compare_rev(i) = M3_rev;
end

%%
xi_e_set = 0.2:0.2:1.6;

figure;
plot(xi_e_set,M0_e_compare_cost,'-x','linewidth',2,'markersize',16);
hold on
plot(xi_e_set,M1_e_compare_cost,'-s','linewidth',2,'markersize',12);
hold on
plot(xi_e_set,M2_e_compare_cost,'-*','linewidth',2,'markersize',16);
hold on
plot(xi_e_set,M3_e_compare_cost,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks(0.2:0.2:1.6);
xlabel('\xi_e','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Total cost','FontSize',18);

%%
Psi_set = 0.2:0.2:1.6;

figure;
plot(xi_e_set,M0_e_compare_rev,'-x','linewidth',2,'markersize',16);
hold on
plot(xi_e_set,M1_e_compare_rev,'-s','linewidth',2,'markersize',12);
hold on
plot(xi_e_set,M2_e_compare_rev,'-*','linewidth',2,'markersize',16);
hold on
plot(xi_e_set,M3_e_compare_rev,'-v','linewidth',2,'markersize',12);
grid on
set(gca,'FontSize',18);

xlim([0.2 1.6]);
xticks(0.2:0.2:1.6);
xlabel('\xi_e','FontSize',18);
legend('M0','M1','M2','M3','FontSize',16)
ylabel('Revenue','FontSize',18);

