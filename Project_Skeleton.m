%% Project Skeleton Code %%

clear; clc; close all;


%% initial parameters

%link lengths, units of metres
r1 = 7.8*10^(-2); % o2o3
r2 = 2.5*10^(-2); % o2a2
r3 = 13.8*10^(-2); % o3B
r5 = 4.75*10^(-2); % BC
r7 = 17.1*10^(-2); % O4O3

syms t theta2(t) theta3(t) r4(t) theta5(t) r6(t); %establishes variables are all functions of t

% link 2 movement
dtheta2 = 1800; % deg/s from 300 rpm
ddtheta2 = 0; %#ok<NASGU> % constant angular velocity
theta2(t) = t*dtheta2 % theta2 value 

%% Part 1- Calculations for kinematic variables, using LCEs
% Found using link closure equations as detailed in the report

theta3(t) = atand((r2*sind(theta2(t)))/(r2*cosd(theta2(t))-r1))+180

r4(t) = (r2*cosd(theta2(t))-r1)/cosd(theta3(t))

theta5(t) = acosd((r7+r3*cosd(theta3(t)))/r5(t))

r6(t) = r5(t)*sind(theta5(t))-r3*sind(theta3(t))

%% Derivative equations of kinematic vars (d/dt) 
syms dtheta3 ddtheta3 dr4 ddr4 dtheta5 ddtheta5 dr6

dtheta3(t) = diff(theta3(t), t); % differentiates theta3_eqn with respect to t
ddtheta3(t) = diff(theta3(t), t, 2); % differentiates dtheta3_eqn with respect to t

dr4(t) = diff(r4(t), t);
ddr4(t) = diff(r4(t), t, 2);

dtheta5(t) = diff(theta5(t), t);
ddtheta5(t) = diff(theta5(t), t, 2);

dr6(t) = diff(r6(t), t);
ddr6(t) = diff(r6(t), t, 2);

%% Calculate Values %%

theta3_array = ones(20);

for s = 0:20 % calculates values over the course of 20 seconds
%     theta2_num = t*dtheta2; % establish numerical value of theta2
    theta3_array(s) = subs(theta3_eqn, t, 1);
%     theta3_num = vpa(180 + solve(subs(theta3_eqn, t, s), theta3))
  
    %theta3_array(deg) = theta3_num
end

%% Plot vars;

% Plot all desired deliverables. 
% 
% figure (1)
% plot(theta2_num, theta3_num)
% grid on;
% title('$\theta_3$ vs $\theta_2$', 'Interpreter','latex')
% xlabel('\theta_2   unit: degree')
% ylabel('\theta_3   unit: degree')
 
% *****************************************************
%% Part 2 - Force and Moment Calculation

%%initial parameters:

dtheta2 = 1800; % theta2 dot
ddtheta2 = 0; % theta2 doble-dot - second derivative

rho = 2.7/1000; % density, kg/cm3
r = 0.25; % radius, cm

m2 = pi*(r^2)*rho*r2*100 ; % link 2, O2A2, kg
m3 = pi*(r^2)*rho*r3*100 ; % link 3, O3B, kg
m5 = pi*(r^2)*rho*r5*100 ; % link 5, BC, kg
m4 = 5/1000 ; % slider 4, kg
m6 = 5/1000 ; % slider 6, kg
%Formula: I = 1/12*m*L^2;
IG5 = m5*(r5^2)/12; %Moment of Inertia, link 5, kg*m^2
IG3 = m3*(r3^2)/12; %MOI, link 3, kg*m^2
IG2 = m2*(r2^2)/12; %MOI, link 2, kg*m^2


M12_list = [];
theta2_list = [];
Fs_list = [];  % shaking force
alpha_s_list = []; % direction of a shaking force
Ms_list =[]; % Shaking moment
Fij_list = []; % Forces
Fij_alpha = []; % Angles at which forces are acting


for theta2 = 0:1:360

    % kinematic variables are caculated based on loop eqn
    r3 = % ENTER YOUR CODE HERE %;
    dr3 = % ENTER YOUR CODE HERE %;
    ddtheta3 = % ENTER YOUR CODE HERE %;

% and so on    

    B = get_ma_vector(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);
    
    A = get_A_matrix(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);

    x = A\ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
    % M12:
    M12 = x(% ENTER YOUR CODE HERE%);
    M12_list = [M12_list; M12];
    
    Fijx = x(% ENTER YOUR CODE HERE%);
    Fijy = x(% ENTER YOUR CODE HERE%);
    
    % Magnitudes of all forces: 
    % Atan is defined on [-pi/2; pi/2]. 
    % This if clause will help to adjust the value of the angle 
    % to its true value:	
    Fij_list = [Fij_list; % ENTER YOUR CODE HERE%];

    
    % Directions of all forces:    
    fx = % ENTER YOUR CODE HERE%;
    fy = % ENTER YOUR CODE HERE%;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    Fij_alpha = [Fij_alpha; alpha_f];

    % and so on
    
  
    % Collecting the values of theta2:
    theta2_list = [theta2_list, theta2];
     
   
    
end


% Regular and Polar plots:
% Might have to transpose the Force vectors for polar plot. Do so if needed
% Polar plot only works with radians so will have to do it accordingly

figure (3)
plot(theta2_list,M12_list)
grid on;
title('M_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')


% Convert degrees to the radians
theta2_rad = deg2rad(theta2_list);

figure (4)
polarplot(Fij_alpha,Fij_list)
grid on;
title('F_{ij} polar plot')

% and so on ...
