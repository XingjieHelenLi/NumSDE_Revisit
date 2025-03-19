%function moment_evolution()
% trajectories of 1st and 2nd moments for each method

close all

%% Model setting


eta =1.4; 
% eta = 0.1;
x0 = 1;


%% numerical setting

T = 20;
h_array = [0.001, 0.01, 0.1];
h_L = length(h_array);

h = h_array(1);  t_array = (0:floor(T/h))*h;
figure(77);
subplot(1,3,1)
semilogy(t_array, x0*exp(-t_array),'k-','DisplayName','Analytic');
hold on;

subplot(1,3,2)
semilogy(t_array, x0*exp(-t_array),'k-','DisplayName','Analytic');
hold on;

subplot(1,3,3)
semilogy(t_array, x0*exp(-t_array),'k-','DisplayName','Analytic');
hold on;



if eta ==1
    ext_Mom2nd = x0^2*exp(-t_array)+ 1-exp(-t_array)+2*x0*t_array.*exp(-t_array);
else
    ext_Mom2nd = x0^2*exp(-(2-eta^2)*t_array)+1/(2-eta^2)*(1-exp(-(2-eta^2)*t_array))...
        +2*eta*x0/(1-eta^2).*(exp(-t_array)-exp(-(2-eta^2)*t_array));
end

figure(88);
subplot(1,3,1)
semilogy(t_array, ext_Mom2nd,'k-','DisplayName','Analytic');
hold on;

subplot(1,3,2)
semilogy(t_array, ext_Mom2nd,'k-','DisplayName','Analytic');
hold on;

subplot(1,3,3)
semilogy(t_array, ext_Mom2nd,'k-','DisplayName','Analytic');
hold on;

color_array =['b','r','g','m','c'];
for h_idx =1:h_L
    h = h_array(h_idx);
    evl_h_ratio = 1-h;
    evl_h_ratio_sq = evl_h_ratio.^2;

    Nt = floor(T/h)+1;
    Mom1st_EM = zeros(Nt+1,1);
    Mom2nd_EM = Mom1st_EM;

    %% EM method
    Mom1st_EM(1) = x0;
    Mom2nd_EM(1) = x0^2;

    for i =1:Nt
        Mom1st_EM(i+1)= evl_h_ratio*Mom1st_EM(i);
        Mom2nd_EM(i+1) = ((1-h)^2+h*eta^2)*Mom2nd_EM(i)+2*h*eta*Mom1st_EM(i)+h;
    end

    switch h_idx
        case 1
            figure(77); subplot(1,3,1); semilogy([0:Nt]*h, Mom1st_EM,'--','color',color_array(h_idx),'DisplayName','EM');
            figure(88); subplot(1,3,1); semilogy([0:Nt]*h, Mom2nd_EM,'--','color',color_array(h_idx),'DisplayName','EM');
        case 2
            figure(77); subplot(1,3,2); semilogy([0:Nt]*h, Mom1st_EM,'--','color',color_array(h_idx),'DisplayName','EM');
            figure(88); subplot(1,3,2); semilogy([0:Nt]*h, Mom2nd_EM,'--','color',color_array(h_idx),'DisplayName','EM');
        case 3
            figure(77); subplot(1,3,3); semilogy([0:Nt]*h, Mom1st_EM,'--','color',color_array(h_idx),'DisplayName','EM');
            figure(88); subplot(1,3,3); semilogy([0:Nt]*h, Mom2nd_EM,'--','color',color_array(h_idx),'DisplayName','EM');
    end
           
    %% Milstein method

    Mom1st_Mil = zeros(Nt+1,1);
    Mom2nd_Mil = Mom1st_Mil;
    Mil_h_ratio = h+0.5*h^2*eta^2;

    Mom1st_Mil(1) = x0;
    Mom2nd_Mil(1) = x0^2;

    for i =1: Nt
        Mom1st_Mil(i+1) = evl_h_ratio*Mom1st_Mil(i);
        Mom2nd_Mil(i+1) = evl_h_ratio_sq*Mom2nd_Mil(i)+(Mil_h_ratio)*(1+2*eta*Mom1st_Mil(i) +eta^2*Mom2nd_Mil(i));
    end

    switch h_idx
        case 1
            figure(77); subplot(1,3,1); semilogy([0:Nt]*h, Mom1st_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
            figure(88); subplot(1,3,1); semilogy([0:Nt]*h, Mom2nd_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
        case 2
            figure(77); subplot(1,3,2); semilogy([0:Nt]*h, Mom1st_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
            figure(88); subplot(1,3,2); semilogy([0:Nt]*h, Mom2nd_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
        case 3
            figure(77); subplot(1,3,3); semilogy([0:Nt]*h, Mom1st_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
            figure(88); subplot(1,3,3); semilogy([0:Nt]*h, Mom2nd_Mil,'-.','color',color_array(h_idx),'DisplayName','Milstein');
    end

    %% Stochastic Heun method

    Mom1st_SH = zeros(Nt+1,1);
    Mom2nd_SH = Mom1st_SH;

    SH_h_ratio1 = 1-h+1/8*h^2*(2+eta^2)^2;
    SH_h_ratio2 = 1/8*eta*h^2*(2+eta^2);

    SH_h_ratio3 = h-1/2*h^2*(2+eta^2)+1/4*h^3*(1+eta^2)^2+1/64*h^3*eta^2*(2+eta^2)^2;
    SH_h_ratio4 = 2*h*eta-1/4*h^2*eta*(10+3*eta^2)+1/4*h^3*eta*(2+5*eta^2+2*eta^4)...
                    +1/32*h^4*eta*(2+eta^2)^3;
    SH_h_ratio5 = 1-h*(2-eta^2)-1/4*h^2*((2+eta^2)^2-12)+1/4*h^3*(eta^6+3*eta^4-4)+1/64*h^4*(2+eta^2)^4;


    Mom1st_SH(1) = x0;
    Mom2nd_SH(1) = x0^2;

    for i =1: Nt
        Mom1st_SH(i+1) = SH_h_ratio1*Mom1st_SH(i)+SH_h_ratio2;        
        Mom2nd_SH(i+1) = SH_h_ratio3+SH_h_ratio4*(Mom1st_SH(i))+SH_h_ratio5*Mom2nd_SH(i);
    end

    switch h_idx
        case 1
            figure(77); subplot(1,3,1); semilogy([0:Nt]*h, Mom1st_SH,':','color',color_array(h_idx),'DisplayName','SH');
            figure(88); subplot(1,3,1); semilogy([0:Nt]*h, Mom2nd_SH,':','color',color_array(h_idx),'DisplayName','SH');
        case 2
            figure(77); subplot(1,3,2); semilogy([0:Nt]*h, Mom1st_SH,':','color',color_array(h_idx),'DisplayName','SH');
            figure(88); subplot(1,3,2); semilogy([0:Nt]*h, Mom2nd_SH,':','color',color_array(h_idx),'DisplayName','SH');
        case 3
            figure(77); subplot(1,3,3); semilogy([0:Nt]*h, Mom1st_SH,':','color',color_array(h_idx),'DisplayName','SH');
            figure(88); subplot(1,3,3); semilogy([0:Nt]*h, Mom2nd_SH,':','color',color_array(h_idx),'DisplayName','SH');
    end
    
    %% 3-stage Runge Kutta (RK3)

    Mom1st_RK3 = zeros(Nt+1,1);
    Mom2nd_RK3 = Mom1st_RK3;
    
    RK3_h_ratio1 = -1/24*h^2*eta*(2+3*eta^2)-1/48*h^3*eta*(2+eta^2)^2;
    RK3_h_ratio2 = 1-h+1/8*h^2*(4-eta^4)-1/48*h^3*(2+eta^2)^3;
    RK3_h_ratio3 = h-1/2*h^2*(2-eta^2)+1/12*h^3*(8-eta^4)-1/192*h^4*(32+20*eta^2-44*eta^4-27*eta^6)...
                   +1/288*h^5*(2+eta^2)^2*(2+7*eta^2+6*eta^4)+1/2304*h^6*eta^2*(2+eta^2)^4;

    RK3_h_ratio4 = 2*h*eta-1/12*h^2*eta*(38-9*eta^2)+1/24*h^3*eta*(56+2*eta^2-5*eta^4)...
                   -1/96*h^4*eta*(72+44*eta^2-50*eta^4-27*eta^6)...
                     +1/72*h^5*eta*(2+eta^2)^3*(1+3*eta^2)+1/1152*h^6*eta*(2+eta^2)^5;

    RK3_h_ratio5 = 1-h*(2-eta^2)+1/4*h^2*(8-8*eta^2+eta^4)-1/24*h^3*(32-36*eta^2+3*eta^6)...
                   +1/192*h^4*(2+eta^2)^2*(28-52*eta^2+27*eta^6)...
                     -1/96*h^5*(2+eta^2)^4*(1-2*eta^2)+1/2304*h^6*(2+eta^2)^6;


    Mom1st_RK3(1) = x0;
    Mom2nd_RK3(1) = x0^2;

    for i =1: Nt
        Mom1st_RK3(i+1) = RK3_h_ratio1+ RK3_h_ratio2*Mom1st_RK3(i);
        Mom2nd_RK3(i+1) = RK3_h_ratio3+ RK3_h_ratio4*Mom1st_RK3(i)+ RK3_h_ratio5*Mom2nd_RK3(i); 

    end

    switch h_idx
        case 1
            figure(77); subplot(1,3,1); semilogy([0:Nt]*h, Mom1st_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
            figure(88); subplot(1,3,1); semilogy([0:Nt]*h, Mom2nd_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
        case 2
            figure(77); subplot(1,3,2); semilogy([0:Nt]*h, Mom1st_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
            figure(88); subplot(1,3,2); semilogy([0:Nt]*h, Mom2nd_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
        case 3
            figure(77); subplot(1,3,3); semilogy([0:Nt]*h, Mom1st_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
            figure(88); subplot(1,3,3); semilogy([0:Nt]*h, Mom2nd_RK3,'.','color',color_array(h_idx),'DisplayName','RK3');
    end



end

title_txt =['eta = ', num2str(eta)];
figure(77);title(title_txt);xlabel('Time'); ylabel('First moment');legend
figure(88);title(title_txt);xlabel('Time'); ylabel('Second moment');legend

 