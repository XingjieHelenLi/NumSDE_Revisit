%function LevelPlot_SDE()
%{
This function is used to compute and plot the stability regions (h vs eta) of first
and second moments for various methods.
%}

close all
%% Parameter settings

N = 1001;


% Step 1: Define the range for x and y
h_array = linspace(0.001, 2.6,N);
eta_array = linspace(0, sqrt(2), N);

% Step 2: Create a grid of points
[X, Y] = meshgrid(eta_array, h_array);

% Step 3: Define the function f(x, y) for 2nd moment stability
F_SH = zeros(N, N);
F_Mil = F_SH;
F_EM = F_SH;
F_IRK3 = F_SH;

F_IRK3_1st = F_IRK3;


%% Second Moment stab constraint

% Euler-Maruyama
EM_stab = @(eta, h) (1-h).^2+h.*eta.^2;

% Milstein 
Mil_stab = @(eta, h) (1-h).^2+h.*eta.^2+ 0.5*h.^2.*eta.^4;

% Stochastic Heun method
const1 = @(eta)  (2-eta.^2);
const2 = @(eta)  (4*eta.^2+eta.^4-8)/4;
const3 = @(eta)  (-4+3*eta.^4+eta.^6)/4;
const4 = @(eta)   (2+eta.^2).^4/64;
SH_stab = @(eta,h)  1-h.*const1(eta)-h.^2.*const2(eta)+h.^3.*const3(eta)+h.^4.*const4(eta);


% IRK3 method

const5 = @(eta) (8-8*eta.^2+eta.^4)/4;
const6 = @(eta) (32-36*eta.^2+3*eta.^6)/24;
const7 = @(eta) (2+eta.^2).^2*(28-52*eta.^2+27*eta.^4)/192;
const8 = @(eta) (2+eta.^2).^4*(1-2*eta.^2)/96;
const9 = @(eta) (2+eta.^2).^6/2304;

IRK3_stab = @(eta, h) 1-h*const1(eta)+h.^2.*const5(eta)-h.^3.*const6(eta) ...
                     +h.^4.*const7(eta)-h.^5.*const8(eta)+h.^6*const9(eta);



for i =1:N
    h = h_array(i);
    for j=1:N

        eta = eta_array(j);
        F_EM(i,j) = EM_stab(eta,h);
        F_Mil(i,j) = Mil_stab(eta, h);
        F_SH(i,j)= SH_stab(eta,h);
        F_IRK3(i,j) = IRK3_stab(eta, h);

    end
end


stab_F_EM = abs(F_EM)<1;
stab_F_Mil = abs(F_Mil)<1;
stab_F_SH = abs(F_SH)<1;
stab_F_IRK3 = abs(F_IRK3)<1;

h_1stMom_cond_EM = 2*ones(size(eta_array));
h_1stMom_cond_Mil = 2*ones(size(eta_array));
h_1stMom_cond_SH = 2./((1+0.5*eta_array.^2)).^2;



%% First moment stability constraint for IRK3
aa = @(eta) (2+eta^2)^3/48;
bb = @(eta) -(4-eta^4)/8;
cc = 1;
dd = -2;

h_1stMom_cond_IRK3 = zeros(size(eta_array));


for i =1:N
    h = h_array(i);
    for j=1:N

        eta = eta_array(j);
        F_IRK3_1st(i,j) = 1-h+1/8*h^2*(4-eta^4)-1/48*h^3*(2+eta^2)^3;

    end
end

figure; contour(X, Y,F_IRK3_1st, [-1 -0.995],'c:','LineWidth',2);
xlabel('eta')
ylabel('h');

axis equal;
title('First moment stability bounary of RK3')

for i = 1:N
    eta = eta_array(i);
    a = aa(eta); b= bb(eta); c= cc; d= dd;
    coefficients = [a b c d];
    roots_of_polynomial = roots(coefficients);
    h_1stMom_cond_IRK3(i) = max(0, max(real(roots_of_polynomial)));
    
end



% Step 4: Plot the contour at the level f(x, y) = -1 and 1
figure;
contourf( X, Y, stab_F_EM);  % Level set where f(x, y) = 1
hold on;
plot(eta_array, h_1stMom_cond_EM,'g-','LineWidth', 2);
xlabel('eta');
ylabel('h');

axis equal;
title('Stability bounary of EM')

figure;
contourf( X, Y, stab_F_Mil);  % Level set where f(x, y) = 1
hold on;
plot(eta_array, h_1stMom_cond_Mil,'m:','LineWidth', 2);
xlabel('eta');
ylabel('h');
%title('Level Set of f(x, y) = x^2 + y^2');
axis equal;
title('Stab boundary of Milstein')

figure;
contourf( X, Y, stab_F_SH);  % Level set where f(x, y) = 1
hold on;
plot(eta_array, h_1stMom_cond_SH,'r--', 'LineWidth', 2)

xlabel('eta');
ylabel('h');
%title('Level Set of f(x, y) = x^2 + y^2');
axis equal;

title('Stability bounary of SH')


figure;
contourf( X, Y, stab_F_IRK3);  % Level set where f(x, y) = 1
hold on;
plot(eta_array, h_1stMom_cond_IRK3,'m:', 'LineWidth', 2)

xlabel('eta');
ylabel('h');
%title('Level Set of f(x, y) = x^2 + y^2');
axis equal;

title('Stability bounary of RK3')

%figure;

%surf(X,Y, F_SH)


figure;





[C1,h1]= contourf( X, Y, stab_F_EM);hold on;

%customColormap1 = [linspace(0.25, 0.75, 256)', linspace(0, 0, 256)', linspace(0.75, 0.25, 256)'];
%colormap(customColormap1);

[C2,h2]= contourf( X, Y, stab_F_Mil);
%customColormap2 = [linspace(0.25, 0.5, 256)', linspace(0.25, 0.25, 256)', linspace(0.5, 0.25, 256)'];
%colormap(customColormap2);
[C3,h3]= contourf( X, Y, stab_F_SH);

[C4,h4]= contourf( X, Y, stab_F_IRK3); 

%customColormap3 = [linspace(0.7, 0.8, 32)', linspace(0.7, 0.8, 32)', linspace(0.7, 0.8, 32)'];
%colormap(customColormap3);

customColormap3 = [linspace(0.5, 0.75, 32)', linspace(0.75, 0.75, 32)', linspace(0.75, 0.5, 32)'];
colormap(customColormap3);


set(h1, 'FaceAlpha', 0.35);
set(h3, 'FaceAlpha', 0.45);
set(h2, 'FaceAlpha', 0.6);
set(h4, 'FaceAlpha', 0.75);

plot(eta_array, h_1stMom_cond_EM,'g-','LineWidth', 2);
plot(eta_array, h_1stMom_cond_Mil,'r:','LineWidth', 2);
plot(eta_array, h_1stMom_cond_SH,'b--', 'LineWidth', 2); 
plot(eta_array, h_1stMom_cond_IRK3,'m:', 'LineWidth', 2)

%colormap("colorcube");  % Use the colormap
%customColormap = [linspace(0, 1, 128)', linspace(0, 0, 128)', linspace(1, 0, 128)'];  % Blue to Red
%colormap(customColormap);
xlabel('eta'); ylabel('h');



%legend('EM','Milstein','SH','RK3')