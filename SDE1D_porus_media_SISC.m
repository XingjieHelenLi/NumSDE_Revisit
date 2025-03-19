%function SDE1D_porus_media_SISC()


close all;clc;

% Load advection velocity drift and noise diffusion coeff
% both velocity v(x) and diffusion sigma(x) are defined on discrete grid stencils x ver domain [0, L]
% we will interpolate the continuous function by first smoothing data and
% then spline interpolation

% load SDE coeff data, v(x) and sigma(x), for LARGE eta 
%load('RandFlow_Small.mat','v','sigma','x','L');


% load SDE coeff data, v(x) and sigma(x), for SMALL eta 
load('RandFlow_Large.mat','v','sigma','x','L');

disp('The modeling parameters for the random Darcy flow are')


% Spline interpolation for sigma(x) and v(x) over domain [0, L]

smooth_sigma = smoothdata(sigma,'gaussian',7);
smooth_v = smoothdata(v,'gaussian',7);

splinefit_sigma = spline(x,smooth_sigma);
firDeriv_spline_sigma = fnder(splinefit_sigma,1);
secDeriv_spline_sigma = fnder(splinefit_sigma,2);

splinefit_v = spline(x,smooth_v);
firDeriv_spline_v = fnder(splinefit_v,1);

figure(77); hold on; %plot(x,sigma);
plot(x,smooth_sigma,'r-','LineWidth',2); %title('diffusion coeff');
ylabel('$\sigma(x)$','FontSize',12,'Interpreter','latex');
xlabel('$x$','FontSize',12,'Interpreter','latex');

figure(88); hold on; %plot(x,v);
plot(x,smooth_v,'b-','LineWidth',2);%title('velocity')

ylabel('$u(x)$','FontSize',12,'Interpreter','latex');
xlabel('$x$','FontSize',12,'Interpreter','latex');

% plot color array
color_array = ['r','b','m','g','y','c'];

% total number of trajectories
tot_traj =1500; 

% simulation time T and dt
T =100;
dt_array =[0.01,2]; 
dtL = length(dt_array);
dt_ratio = [1,200];

dt_fine = dt_array(1);
Nt_fine = floor(T/dt_fine);


% initial position of 1D SDE 
X0 = 0.55*L;

% stochastic terms for fine dt and coarse dt
dWdt_fine = randn(tot_traj, Nt_fine+1)*sqrt(dt_fine);
dWdt_tilde_fine = randn(tot_traj,Nt_fine+1)*sqrt(dt_fine);


% start 1D SDE simulation via different time integrators

for dt_i = 1:dtL

    dt = dt_array(dt_i)
    ratio = dt_ratio(dt_i);


    Nt = floor(T/dt);
    X_EM_ave = zeros(1,Nt+1);
    X_Mil_ave = X_EM_ave;
    X_SH_ave = X_EM_ave;
    X_IRK3_ave = X_EM_ave;

    EM_tot_traj = tot_traj;
    Mil_tot_traj = tot_traj;
    SH_tot_traj = tot_traj;
    IRK3_tot_traj = tot_traj;

    for traj_idx =1: tot_traj
        
        traj_idx
        
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

            % check if the random particle is still in the domain [0,L]
            % if not, then discard this trajectory
            if X_EM_curr<0 || X_EM_curr>=L
                disp(ti);disp(traj_idx);
                disp('EM out of domain');
                EM_tot_traj= EM_tot_traj-1;
                X_EM = 0*X_EM;
                break;

            end

            X_EM(ti+1) = X_EM_curr + ppval(splinefit_v,X_EM_curr).*dt + ppval(splinefit_sigma,X_EM_curr).*dWdt_curr ;

        end

        X_EM_ave = X_EM_ave+X_EM;

        %% Milstein method
        for ti = 1:Nt
            X_Mil_curr   = X_Mil(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));

            % check if the random particle is still in the domain [0,L]
            % if not, then discard this trajectory
            if X_Mil_curr<0 || X_Mil_curr>=L
                disp(ti);disp(traj_idx);
                disp('Milstein out of domain');
                X_Mil = 0*X_Mil;
                Mil_tot_traj= Mil_tot_traj-1;
                break;
            end
            X_Mil(:,ti+1) = X_Mil_curr + ppval(splinefit_v,X_Mil_curr).*dt + ppval(splinefit_sigma,X_Mil_curr).*dWdt_curr ...
                +0.5*ppval(firDeriv_spline_sigma,X_Mil_curr).*(dWdt_curr.^2-dt);
        end

        X_Mil_ave = X_Mil_ave+X_Mil;


        %% Stochastic Heun method
        for ti = 1:Nt

            X_SH_curr   = X_SH(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));

            % check if the random particle is still in the domain [0,L]
            % if not, then discard this trajectory 
            if X_SH_curr<0 || X_SH_curr>=L
                disp(ti);disp(traj_idx);
                disp('SH out of domain');
                SH_tot_traj= SH_tot_traj-1;
                X_SH = 0*X_SH;
                break;
            end
            % stage 1
            F1 = ppval(splinefit_v, X_SH_curr)-0.5*ppval(splinefit_sigma,X_SH_curr)*ppval(firDeriv_spline_sigma,X_SH_curr);
            g1 = ppval(splinefit_sigma,X_SH_curr);
            % stage 2
            X_SH_temp = X_SH_curr+F1*dt+g1.*dWdt_curr;
            F2 = ppval(splinefit_v, X_SH_temp)-0.5*ppval(splinefit_sigma,X_SH_temp)*ppval(firDeriv_spline_sigma,X_SH_temp);
            g2 = ppval(splinefit_sigma,X_SH_temp);
            X_SH(:,ti+1) = X_SH_curr + 0.5*(F1+F2)*dt+0.5*(g1+g2).*dWdt_curr;

        end

        X_SH_ave = X_SH_ave+X_SH;


        %% Improved 3-stage Stochastic Runge-Kutta 3
        for ti = 1:Nt


            X_IRK3_curr   = X_IRK3(ti);
            dWdt_curr = sum(dWdt_fine_traj((ti-1)*ratio+1:ti*ratio));
            dWdt_curr_tilde = sum(dWdt_tilde_fine_traj((ti-1)*ratio+1:ti*ratio));

            % check if the random particle is still in the domain [0,L]
            % if not, then discard this trajectory

            if X_IRK3_curr<0 || X_IRK3_curr>=L
                disp(ti);disp(traj_idx);
                disp('IR3 out of domain');
                IRK3_tot_traj= IRK3_tot_traj-1;
                X_IRK3 = 0*X_IRK3;
                break;
            end

            % Stage 1-3
            F1 = ppval(splinefit_v, X_IRK3_curr)-0.5*ppval(splinefit_sigma,X_IRK3_curr)*ppval(firDeriv_spline_sigma,X_IRK3_curr);
            g1 = ppval(splinefit_sigma,X_IRK3_curr);

            X_IRK3_temp1 = X_IRK3_curr+1/3*F1*dt+1/3*g1.*dWdt_curr;

            F2 = ppval(splinefit_v, X_IRK3_temp1)-0.5*ppval(splinefit_sigma,X_IRK3_temp1)*ppval(firDeriv_spline_sigma,X_IRK3_temp1);
            g2 = ppval(splinefit_sigma,X_IRK3_temp1);

            X_IRK3_temp2 = X_IRK3_curr+2/3*F2*dt+2/3*g2.*dWdt_curr;

            F3 = ppval(splinefit_v, X_IRK3_temp2)-0.5*ppval(splinefit_sigma,X_IRK3_temp2)*ppval(firDeriv_spline_sigma,X_IRK3_temp2);
            g3 = ppval(splinefit_sigma,X_IRK3_temp2);


            X_IRK3_next = X_IRK3_curr +1/4*(F1+3*F3)*dt+1/4*(g1+3*g3).*dWdt_curr ...
                + 1/(2*sqrt(3))*(ppval(firDeriv_spline_v, X_IRK3_curr).*g1 ...
                -ppval(splinefit_v,X_IRK3_curr).*ppval(firDeriv_spline_sigma,X_IRK3_curr)...
                -0.5*ppval(secDeriv_spline_sigma,X_IRK3_curr)*g1.^2)*dt.* dWdt_curr_tilde;
            X_IRK3(ti+1) = X_IRK3_next;

        end

        X_IRK3_ave = X_IRK3_ave+X_IRK3;



    end

    % compute the mean trajectory for each numerical integrator
    if EM_tot_traj>0;    X_EM_ave = X_EM_ave/EM_tot_traj; end
    if Mil_tot_traj>0; X_Mil_ave = X_Mil_ave/Mil_tot_traj; end
    if SH_tot_traj>0; X_SH_ave = X_SH_ave/SH_tot_traj; end
    if IRK3_tot_traj>0; X_IRK3_ave = X_IRK3_ave/IRK3_tot_traj; end

    % record the mean trajectory for fine dt
    if ratio == 1
        X_EM_ave_fine = X_EM_ave;
        X_Mil_ave_fine = X_Mil_ave;
        X_SH_ave_fine= X_SH_ave;
        X_IRK3_ave_fine = X_IRK3_ave;
    end

    % plot the mean trajectory for each dt
    figure(1); hold on;
    plot((0:Nt)*dt, X_EM_ave,'-','Color',color_array(dt_i));
    plot((0:Nt)*dt, X_Mil_ave,':','LineWidth',1.5,'Color',color_array(dt_i));
    plot((0:Nt)*dt, X_SH_ave,'-.','Color',color_array(dt_i));
    plot((0:Nt)*dt, X_IRK3_ave,'--','Color',color_array(dt_i));
    pause(1)
    
    % record ALL fine trajectories
    if ratio==1
        EM_tot_traj_fine = EM_tot_traj; Mil_tot_traj_fine = Mil_tot_traj;
        SH_tot_traj_fine = SH_tot_traj; IRK3_tot_traj_fine=  IRK3_tot_traj;
    end

end

legend('EM','Milstein','SH','IRK3')
xlabel('time');ylabel('domain')


temp_ratio = 1;
figure; hold on
dt= dt_array(1); Nt_fine = floor(T/dt);
plot((0:Nt_fine)*dt, X_EM_ave_fine,'-','Color',color_array(1));
plot((0:Nt_fine)*dt, X_Mil_ave_fine,':','LineWidth',1.5,'Color',color_array(1));
plot((0:Nt_fine)*dt, X_SH_ave_fine,'-.','Color',color_array(1));
plot((0:Nt_fine)*dt, X_IRK3_ave_fine,'--','Color',color_array(1));
dt = dt_array(2);
plot((0:Nt)*dt, X_EM_ave*temp_ratio,'-','Color',color_array(dt_i));
plot((0:Nt)*dt, X_Mil_ave*temp_ratio,':','LineWidth',1.5,'Color',color_array(dt_i));
plot((0:Nt)*dt, X_SH_ave*temp_ratio,'-.','Color',color_array(dt_i));
plot((0:Nt)*dt, X_IRK3_ave*temp_ratio,'--','Color',color_array(dt_i));


figure; hold on
dt= dt_array(1); Nt_fine = floor(T/dt/2);
plot((0:Nt_fine)*dt, X_EM_ave_fine(1:Nt_fine+1),'-','Color',color_array(1),'DisplayName','Euler-Maruyama');
plot((0:Nt_fine)*dt, X_Mil_ave_fine(1:Nt_fine+1),':','LineWidth',1.5,'Color',color_array(1),'DisplayName','Milstein');
plot((0:Nt_fine)*dt, X_SH_ave_fine(1:Nt_fine+1),'-.','Color',color_array(1),'DisplayName','Stochastic Heun');
plot((0:Nt_fine)*dt, X_IRK3_ave_fine(1:Nt_fine+1),'--','Color',color_array(1),'DisplayName','3-stage Runge-Kutta');
dt = dt_array(2);Nt=floor(T/dt/2);
plot((0:Nt)*dt, X_EM_ave(1:Nt+1)*temp_ratio,'-','Color',color_array(dt_i));
plot((0:Nt)*dt, X_Mil_ave(1:Nt+1)*temp_ratio,':','LineWidth',1.5,'Color',color_array(dt_i));
plot((0:Nt)*dt, X_SH_ave(1:Nt+1)*temp_ratio,'-.','Color',color_array(dt_i));
plot((0:Nt)*dt, X_IRK3_ave(1:Nt+1)*temp_ratio,'--','Color',color_array(dt_i));
xlabel('time','FontSize',12);ylabel('position','FontSize',12);
legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta')


fine_ave = (X_EM_ave_fine+X_Mil_ave_fine+X_SH_ave_fine+X_IRK3_ave_fine)/4;
EM_coarse_err = fine_ave(1:ratio:end)-X_EM_ave;
Mil_coarse_err = fine_ave(1:ratio:end)-X_Mil_ave;
SH_coarse_err = fine_ave(1:ratio:end)-X_SH_ave;
IRK3_coarse_err= fine_ave(1:ratio:end)-X_IRK3_ave;

figure;hold on;
Nt=floor(T/dt/1);
plot((0:Nt)*dt, EM_coarse_err(1:Nt+1)*temp_ratio,'-','Color',color_array(1));
plot((0:Nt)*dt, Mil_coarse_err(1:Nt+1)*temp_ratio,':','LineWidth',1.5,'Color',color_array(2));
plot((0:Nt)*dt, SH_coarse_err(1:Nt+1)*temp_ratio,'-.','Color',color_array(3));
plot((0:Nt)*dt, IRK3_coarse_err(1:Nt+1)*temp_ratio,'--','Color',color_array(4));
xlabel('time','FontSize',12);ylabel('position','FontSize',12);
legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta')
xlabel('time','FontSize',12);ylabel('fine - coarse','FontSize',12);


figure;hold on;
Nt=floor(T/dt/1);
plot((0:Nt)*dt, abs(EM_coarse_err(1:Nt+1))*temp_ratio,'-','Color',color_array(1));
plot((0:Nt)*dt, abs(Mil_coarse_err(1:Nt+1))*temp_ratio,':','LineWidth',1.5,'Color',color_array(2));
plot((0:Nt)*dt, abs(SH_coarse_err(1:Nt+1))*temp_ratio,'-.','Color',color_array(3));
plot((0:Nt)*dt, abs(IRK3_coarse_err(1:Nt+1))*temp_ratio,'--','Color',color_array(4));
xlabel('time','FontSize',12);ylabel('position','FontSize',12);
legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta')
xlabel('time','FontSize',12);ylabel('fine - coarse','FontSize',12);

%%%--------------------



X_SH_ave2 = X_SH_ave; 

EM_coarse_err2 = EM_coarse_err;
Mil_coarse_err2 = Mil_coarse_err;
SH_coarse_err2 = fine_ave(1:ratio:end)-X_SH_ave2;
IRK3_coarse_err2 = IRK3_coarse_err;


figure;
dt= dt_array(1); Nt_fine = floor(T/dt/1);
plot((0:Nt_fine)*dt, X_EM_ave_fine(1:Nt_fine+1),'-','Color',color_array(1),'DisplayName','Euler-Maruyama'); hold on;
plot((0:Nt_fine)*dt, X_Mil_ave_fine(1:Nt_fine+1),':','LineWidth',1.5,'Color',color_array(1),'DisplayName','Milstein');
plot((0:Nt_fine)*dt, X_SH_ave_fine(1:Nt_fine+1),'-.','Color',color_array(1),'DisplayName','Stochastic Heun');
plot((0:Nt_fine)*dt, X_IRK3_ave_fine(1:Nt_fine+1),'--','Color',color_array(1),'DisplayName','3-stage Runge-Kutta');
dt = dt_array(dt_i);Nt=floor(T/dt/1);
plot((0:Nt)*dt, X_EM_ave(1:Nt+1)*temp_ratio,'-','Color',color_array(dt_i));
plot((0:Nt)*dt, X_Mil_ave(1:Nt+1)*temp_ratio,':','LineWidth',1.5,'Color',color_array(dt_i));
plot((0:Nt)*dt, X_SH_ave2(1:Nt+1)*temp_ratio,'-.','Color',color_array(dt_i));
plot((0:Nt)*dt, X_IRK3_ave(1:Nt+1)*temp_ratio,'--','Color',color_array(dt_i));
xlabel('time','FontSize',12);ylabel('position','FontSize',12);
legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta')



figure;hold on;
Nt=floor(T/dt/1);
plot((0:Nt)*dt, abs(EM_coarse_err2(1:Nt+1))*temp_ratio,'-','Color',color_array(1));
plot((0:Nt)*dt, abs(Mil_coarse_err2(1:Nt+1))*temp_ratio,':','LineWidth',1.5,'Color',color_array(2));
plot((0:Nt)*dt, abs(SH_coarse_err2(1:Nt+1))*temp_ratio,'-.','Color',color_array(3));
plot((0:Nt)*dt, abs(IRK3_coarse_err2(1:Nt+1))*temp_ratio,'--','Color',color_array(4));
xlabel('time','FontSize',12);ylabel('position','FontSize',12);
legend('Euler Maruyama','Milstein','Stochastic Heun','3-stage Runge-Kutta')
xlabel('time','FontSize',12);ylabel('|fine - coarse|','FontSize',12);