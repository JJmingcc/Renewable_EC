function [rev, cost,P_G,P_S,c,x,q,uti_percent] = Model_3(I,J,T,alpha,lambda_fore,phi,C_max,rho,e,a,P_idle,P_peak,P_R,P_G_max,delta,theta_grid,Delta_T,gamma_max,E_usage)
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
    variable P_S(J,T)
    variable P_U(J,T)
    variable EM(J,T)
    % master obj funtion
    variable C_u
    variable C_e
    variable C_c
    variable obj
    minimize obj
    subject to 
          obj == C_u + C_c + C_e;
%             obj == C_u + C_c;
            % Demand supply balance: Eq(2)
            for t = 1:T
                alpha * sum(x(:,:,t),2) + q(:,t) >= lambda_fore(:,t); 
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
            for t = 1:T
                for j = 1:J
                    P_G(j,t) + P_R(j,t) >= P_U(j,t) + P_S(j,t);
                    P_G(j,t) <= P_G_max(j,t);
                end
            end

            c >= 0;
            x >= 0;
            q >= 0;
            P_G >= 0;
            P_S >= 0;
cvx_end
cost = full(obj);
P_G = full(P_G);
P_S = full(P_S);
c = full(c);
x = full(x);
q = full(q);
uti_percent = full(workload)/(rho * c);
rev = sum(sum(a.* P_S));

