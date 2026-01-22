function SDE_nonlinear_erf()

%{
This function aims to simulate the following 1D SDE by four numerical
integrators

SDE: dx = f(x)dt+ g(x)dW_t for -L/2 <= x<= L/2
     with reflection boundary condition

f(x) = E*erf((x-x_star2)/D)
g(x) = A*erf((x-x_star)/B)+C

Four schemes:
EM: Euler-Maruyama
Mil: Milstein
SH: Stochastic Heun 
RK3: Improved 3-stage Stochastic Runge Kutta
%}

close all; clc; clear all

% length of physical domain 
L = 16.0;

% eta_case =1: small Eta; eta_case = 2: large eta

eta_case = 2
switch eta_case
    case 1
% small Eta case
    E = -2.0 ;    D = 1.0;     A = 0.3;     B = 4.0;     x_star = 1.0;    C = 1.0; x_star2 = x_star;  
    case 2
% large Eta case
    E = -2.0 ;    D = 1.0;     A = 1.5;     B = 1.0;     x_star = 1.0;    C = 1.0; x_star2 = x_star;
end

%% Plot the drift and diffusion coefficients of the underlying nonlinear pblem
%plot_porus_coeff(L,A,B,C,D,E,x_star, x_star2) 

%% drift f(x), diffusion g(x) and their derivatives

% Define functions
f_x = @(x) E * erf((x-x_star2) / D);
g_x = @(x) A * erf((x - x_star) / B) + C;

% Compute derivatives
der_f_x = @(x) E / D * 2 / sqrt(pi) * exp(-((x-x_star2) / D).^2);
der_g_x = @(x) A / B * 2 / sqrt(pi) * exp(-((x - x_star) / B).^2);
der_2nd_g_x =  @ (x) A/(B^3) * 2/sqrt(pi) * exp(-((x - x_star).^2) ./ (B^2)) .* (-2*(x - x_star));



% compute the exact stationary pdf, first and second moments, and variance
% by the trapzoid rule 

[xx2, rho, mean_x, second_moment, var_x, ~] = ...
    stationary_density_reflecting(L, E, D, A, B, x_star, C, x_star2);

% total number of sample paths
tot_traj =40000; 
% total simulation time for one trajectory
T = 20; 

dt_ratio = 2.^[0,1,2,3,4,5,6]; 
dt_array =0.01*dt_ratio; 
dtL = length(dt_array);

color_array = ['r','b','m','g','y','c','k','r'];
dt_fine = dt_array(1);
Nt_fine = floor(T/dt_fine);

% generate random noise terms for all trajectories for finest time step
% size
dWdt_fine = randn(tot_traj, Nt_fine+1)*sqrt(dt_fine);
dWdt_tilde_fine = randn(tot_traj,Nt_fine+1)*sqrt(dt_fine);

% variables for saving first moment
X_EM_long_mean = zeros(1,dtL);
X_Mil_long_mean = zeros(1,dtL);
X_SH_long_mean = zeros(1,dtL);
X_IRK3_long_mean = zeros(1,dtL);

% variables for saving second moment
X_EM_long_var = zeros(1,dtL);   
X_Mil_long_var = zeros(1,dtL);  
X_SH_long_var = zeros(1,dtL);  
X_IRK3_long_var = zeros(1,dtL); 


%% Boundary values
a= -L/2; b = L/2;

for dt_i = 1:dtL

    dt = dt_array(dt_i)
    ratio = dt_ratio(dt_i);

    % given time step size dt, initialize the varialbes to store the
    % trajectories of each scheme given dt
    Nt = floor(T/dt);
    X_EM_ave = zeros(1,Nt+1);
    X_Mil_ave = X_EM_ave;
    X_SH_ave = X_EM_ave;
    X_IRK3_ave = X_EM_ave;

    EM_tot_traj = tot_traj;
    Mil_tot_traj = tot_traj;
    SH_tot_traj = tot_traj;
    IRK3_tot_traj = tot_traj;

    var_len = max(floor(0.1/dt),5);

    for traj_idx =1: tot_traj
        traj_idx

       
        X0 = mean_x+1.5; % sample the initial data from the stationary mean 

        %dWdt = dW_fine(traj_idx,1:ratio:Nt_fine+1)*sqrt(dt);
        %dWdt_tilde = dW_tilde_fine(traj_idx, 1:ratio:Nt_fine+1)*sqrt(dt);
        dWdt_fine_traj = dWdt_fine(traj_idx,:);
        dWdt_tilde_fine_traj = dWdt_tilde_fine(traj_idx,:);

        X_EM = zeros(1,Nt+1);
        X_EM(1) = X0;

        X_Mil = X_EM;
        X_SH = X_EM;
        X_IRK3 = X_EM;

        %% Euler Maruyama method

        for ti = 1:Nt
            X_EM_curr = X_EM(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));

            X_next = X_EM_curr + f_x(X_EM_curr).*dt + g_x(X_EM_curr).*dWdt_curr ;

            % Reflecting boundary condition
            if X_next < a
                % reflect: position below a => reflect into [a,b]
                X_next = a + (a - X_next);
                if X_next > b
                    % if overshoot double reflect
                    X_next = b - (X_next - b);
                end
            elseif X_next > b
                % reflect from upper boundary
                X_next = b - (X_next - b);
                if X_next < a
                    % if overshoot double reflect
                    X_next = a + (a - X_next);
                end
            end

            X_EM(ti+1) = X_next;

        end

        X_EM_ave = X_EM_ave+X_EM;

        %X_EM_long_mom2(dt_i) = X_EM_long_mom2(dt_i) +X_EM(end)/tot_traj;

        X_EM_long_var(dt_i) = X_EM_long_var(dt_i) + var(X_EM(end-var_len:end))/tot_traj;

        % Milstein method
        for ti = 1:Nt
            X_Mil_curr   = X_Mil(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));

            X_next = X_Mil_curr + f_x(X_Mil_curr).*dt + g_x(X_Mil_curr).*dWdt_curr ...
                +0.5* g_x(X_Mil_curr)*der_g_x(X_Mil_curr).*(dWdt_curr.^2-dt);

            % Reflecting boundary condition
            if X_next < a
                % reflect: position below a => reflect into [a,b]
                X_next = a + (a - X_next);
                if X_next > b
                    % if overshoot double reflect
                    X_next = b - (X_next - b);
                end
            elseif X_next > b
                % reflect from upper boundary
                X_next = b - (X_next - b);
                if X_next < a
                    % if overshoot double reflect
                    X_next = a + (a - X_next);
                end
            end

            X_Mil(ti+1) = X_next;

        end

        X_Mil_ave = X_Mil_ave+X_Mil;

               
         X_Mil_long_var(dt_i) = X_Mil_long_var(dt_i) + var(X_Mil(end-var_len:end))/tot_traj;

        % Stochastic Heun method
        for ti = 1:Nt

            % Stochastic Heun method
            X_SH_curr   = X_SH(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));

            % stage 1
            F1 = f_x(X_SH_curr)-0.5*g_x(X_SH_curr)*der_g_x(X_SH_curr);
            g1 = g_x(X_SH_curr);
            % stage 2
            X_SH_temp = X_SH_curr+F1*dt+g1.*dWdt_curr;
            F2 = f_x(X_SH_temp)-0.5*g_x(X_SH_temp)*der_g_x(X_SH_temp);
            g2 = g_x(X_SH_temp);
            X_next = X_SH_curr + 0.5*(F1+F2)*dt+0.5*(g1+g2).*dWdt_curr;


            % Reflecting boundary condition
            if X_next < a
                % reflect: position below a => reflect into [a,b]
                X_next = a + (a - X_next);
                if X_next > b
                    % if overshoot double reflect
                    X_next = b - (X_next - b);
                end
            elseif X_next > b
                % reflect from upper boundary
                X_next = b - (X_next - b);
                if X_next < a
                    % if overshoot double reflect
                    X_next = a + (a - X_next);
                end
            end

            X_SH(:,ti+1) = X_next;

        end

        X_SH_ave = X_SH_ave+X_SH;

        X_SH_long_var(dt_i) = X_SH_long_var(dt_i) + var(X_SH(end-var_len:end))/tot_traj;

        % Improved 3-stage Stochastic Runge-Kutta 3
        for ti = 1:Nt


            X_IRK3_curr   = X_IRK3(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));
            dWdt_curr_tilde = sum(dWdt_tilde_fine_traj((ti-1)*ratio+1:ti*ratio));

            % Stage 1-3
            F1 = f_x(X_IRK3_curr)-0.5*g_x(X_IRK3_curr)*der_g_x(X_IRK3_curr);
            g1 = g_x(X_IRK3_curr);

            X_IRK3_temp1 = X_IRK3_curr+1/3*F1*dt+1/3*g1.*dWdt_curr;

            F2 = f_x(X_IRK3_temp1)-0.5*g_x(X_IRK3_temp1)*der_g_x(X_IRK3_temp1);
            g2 =g_x(X_IRK3_temp1);

            X_IRK3_temp2 = X_IRK3_curr+2/3*F2*dt+2/3*g2.*dWdt_curr;

            F3 = f_x(X_IRK3_temp2)-0.5*g_x(X_IRK3_temp2)*der_g_x(X_IRK3_temp2);
            g3 = g_x(X_IRK3_temp2);


            X_next = X_IRK3_curr +1/4*(F1+3*F3)*dt+1/4*(g1+3*g3).*dWdt_curr ...
                + 1/(2*sqrt(3))*(der_f_x( X_IRK3_curr).*g1 ...
                -f_x(X_IRK3_curr).*der_g_x(X_IRK3_curr)...
                -0.5*der_2nd_g_x(X_IRK3_curr)*g1.^2)*dt.* dWdt_curr_tilde;

            % Reflecting boundary condition
            if X_next < a
                % reflect: position below a => reflect into [a,b]
                X_next = a + (a - X_next);
                if X_next > b
                    % if overshoot double reflect
                    X_next = b - (X_next - b);
                end
            elseif X_next > b
                % reflect from upper boundary
                X_next = b - (X_next - b);
                if X_next < a
                    % if overshoot double reflect
                    X_next = a + (a - X_next);
                end
            end


            X_IRK3(ti+1) = X_next;

        end

        X_IRK3_ave = X_IRK3_ave+X_IRK3;

        X_IRK3_long_var(dt_i) = X_IRK3_long_var(dt_i) + var(X_IRK3(end-var_len:end))/tot_traj;

    end


    if EM_tot_traj>0;    X_EM_ave = X_EM_ave/EM_tot_traj; end
    if Mil_tot_traj>0; X_Mil_ave = X_Mil_ave/Mil_tot_traj; end
    if SH_tot_traj>0; X_SH_ave = X_SH_ave/SH_tot_traj; end
    if IRK3_tot_traj>0; X_IRK3_ave = X_IRK3_ave/IRK3_tot_traj; end

    if ratio == 1
        X_EM_ave_fine = X_EM_ave;
        X_Mil_ave_fine = X_Mil_ave;
        X_SH_ave_fine= X_SH_ave;
        X_IRK3_ave_fine = X_IRK3_ave;
    end

    figure(111); hold on;
    plot((0:Nt)*dt, X_EM_ave,'-','Color',color_array(dt_i));
    plot((0:Nt)*dt, X_Mil_ave,':','LineWidth',1.5,'Color',color_array(dt_i));
    plot((0:Nt)*dt, X_SH_ave,'-.','Color',color_array(dt_i));
    plot((0:Nt)*dt, X_IRK3_ave,'--','Color',color_array(dt_i));
    pause(1)

    if ratio==1
        EM_tot_traj_fine= EM_tot_traj; Mil_tot_traj_fine= Mil_tot_traj;
        SH_tot_traj_fine= SH_tot_traj; IRK3_tot_traj_fine=  IRK3_tot_traj;
    end
    mean_len = max(floor(0.1/dt),5);
    X_EM_long_mean(dt_i) = mean(X_EM_ave(end-mean_len:end));
    X_Mil_long_mean(dt_i) = mean(X_Mil_ave(end-mean_len:end));
    X_SH_long_mean(dt_i) = mean(X_SH_ave(end-mean_len:end));
    X_IRK3_long_mean(dt_i) = mean(X_IRK3_ave(end-mean_len:end));



end

legend('EM','Milstein','SH','IRK3')
xlabel('time');ylabel('domain')


temp_ratio = 1; 


fine_ave = (X_EM_ave_fine+X_Mil_ave_fine+X_SH_ave_fine+X_IRK3_ave_fine)/4;
EM_coarse_err = fine_ave(1:ratio:end)-X_EM_ave;
Mil_coarse_err = fine_ave(1:ratio:end)-X_Mil_ave;
SH_coarse_err = fine_ave(1:ratio:end)-X_SH_ave;
IRK3_coarse_err= fine_ave(1:ratio:end)-X_IRK3_ave;

dt = dt_array(dtL);

% 
% %%%--------------------



X_SH_ave2 = X_SH_ave; 
EM_coarse_err2 = EM_coarse_err;
Mil_coarse_err2 = Mil_coarse_err;
SH_coarse_err2 = SH_coarse_err;
IRK3_coarse_err2 = IRK3_coarse_err;


figure;hold on;
Nt=floor(T/dt/1);
plot((0:Nt)*dt, abs(EM_coarse_err2(1:Nt+1))*temp_ratio,'-','Color',color_array(1));
plot((0:Nt)*dt, abs(Mil_coarse_err2(1:Nt+1))*temp_ratio,':','LineWidth',1.5,'Color',color_array(2));
plot((0:Nt)*dt, abs(SH_coarse_err2(1:Nt+1))*temp_ratio,'-.','Color',color_array(3));
plot((0:Nt)*dt, abs(IRK3_coarse_err2(1:Nt+1))*temp_ratio,'--','Color',color_array(4));

lgd = legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta');
xlabel('Time','FontSize',12);ylabel('Abs diff in first moment','FontSize',14);
set(lgd,'FontSize',14)

disp('EM long time mean'); disp(X_EM_long_mean);
disp('Mil long time mean'); disp(X_Mil_long_mean);
disp('SH long time mean'); disp(X_SH_long_mean);
disp('IRK3 long time mean'); disp(X_IRK3_long_mean);

disp('EM long time var'); disp(X_EM_long_var);
disp('Mil long time var'); disp(X_Mil_long_var);
disp('SH long time var'); disp(X_SH_long_var);
disp('IRK3 long time var'); disp(X_IRK3_long_var)

exact_long_mean =  mean_x;
figure;
plot_idx = 1;
loglog(dt_array(dtL:-1:plot_idx), abs(X_EM_long_mean(dtL:-1:plot_idx)-exact_long_mean ),'o-','linewidth',2,'color','r','markersize',10);
hold on
loglog(dt_array(dtL:-1:plot_idx), abs(X_Mil_long_mean(dtL:-1:plot_idx)-exact_long_mean ),'s-.','linewidth',2,'color','b','markersize',10);
loglog(dt_array(dtL:-1:plot_idx), abs(X_SH_long_mean(dtL:-1:plot_idx)-exact_long_mean ),'v:','linewidth',2,'color',[1,0.5,0],'markersize',10);
loglog(dt_array(dtL:-1:plot_idx), abs(X_IRK3_long_mean(dtL:-1:plot_idx)-exact_long_mean),'^-.','linewidth',2,'color','g','markersize',10);
% plot ref slope = 1
loglog(dt_array(dtL:-1:plot_idx), (dt_array(dtL:-1:plot_idx)).^(1)/20, '--','linewidth',2,'color','k')

% plot ref slope = 0.5
%loglog(dt_array(plot_idx+3:-1:plot_idx), (dt_array(plot_idx+3:-1:plot_idx)).^(0.5)/5, ':','linewidth',2,'color','k')

xlabel('$h$','interpret','latex','FontSize',14); ylabel('Error of $\mu_T^{(1)}$','interpret','latex','FontSize',14);
lgd = legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta','Ref slope =1');
%lgd = legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta','Ref slope =1','Ref slope = 0.5');

set(lgd,'FontSize',14)



% % Plot the error of numerical values verse exact variance for each scheme
% 
% exact_long_var = var_x;
% 
% figure;
% plot_idx = 1;
% 
% loglog(dt_array(dtL:-1:plot_idx), abs(X_EM_long_var(dtL:-1:plot_idx)-exact_long_var ),'o-','linewidth',2,'color','r','markersize',10);
% hold on
% loglog(dt_array(dtL:-1:plot_idx), abs(X_Mil_long_var(dtL:-1:plot_idx)-exact_long_var ),'s-.','linewidth',2,'color','b','markersize',10);
% loglog(dt_array(dtL:-1:plot_idx), abs(X_SH_long_var(dtL:-1:plot_idx)-exact_long_var ),'v:','linewidth',2,'color',[1,0.5,0],'markersize',10);
% loglog(dt_array(dtL:-1:plot_idx), abs(X_IRK3_long_var(dtL:-1:plot_idx)-exact_long_var),'^-.','linewidth',2,'color','g','markersize',10);
% loglog(dt_array(dtL:-1:plot_idx), (dt_array(dtL:-1:plot_idx)).^(1)/5, '--','linewidth',2,'color','k')
% 
% xlabel('$h$','interpret','latex','FontSize',14); ylabel('Error of $\mathrm{Var}_{T}$','interpret','latex','FontSize',14);
% lgd = legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta','Ref slope =1');
% set(lgd,'FontSize',14)
