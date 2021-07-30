%% Project Skeleton Code %%

clear; clc; close all;

%% initial parameters

%link lengths, units of metres
r1 = 7.8*10^(-2); % o2o3
r2 = 2.5*10^(-2); % o2a2
r3 = 13.8*10^(-2); % o3B
r5 = 4.75*10^(-2); % BC
r7 = 17.1*10^(-2); % o4o3

syms t theta2(t) theta3(t) r4(t) theta5(t) r6(t); %establishes variables are all functions of t

% link 2 movement
dtheta2 = 1800; % deg/s from 300 rpm
ddtheta2 = 0; % constant angular velocity
theta2(t) = t*dtheta2 % theta2 value 

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, using LCEs
% Found using link closure equations as detailed in the report
theta3_eqn = theta3(t) == atand((-sind(theta2(t)))/(3.12 + cosd(theta2(t))))

r4_eqn = r4 == (2.5*sind(theta2)) / (sind(theta3));

theta5_eqn = theta5 == acosd((r7 - r3*cosd(180-theta3)) / (r5)) + 180;
r6_eqn = r6 == r3*sind(180-theta3) - r5*sind(theta5-180);

% Hint: Check if the angle needs to be adjusted to its true value
% Hint: Check this for all other angles too

%% Derivative equations of kinematic vars (d/dt) 
syms dtheta3 ddtheta3 dr4 dtheta5 ddtheta5 dr6

dtheta3_eqn = diff(theta3_eqn, t); % differentiates theta3_eqn with respect to t
ddtheta3_eqn = diff(dtheta3_eqn, t); % differentiates dtheta3_eqn with respect to t
% 
% dr4_eqn = diff(r4_eqn); % dr4(theta2) represented as diff(r4(theta2), theta2)
% 
% dtheta5_eqn = diff(theta5_eqn); % dtheta5(theta2) represented as diff(theta5(theta2), theta2))
% ddtheta5_eqn = diff(dtheta5_eqn); % ddtheta5(theta2) represented as diff(theta5(theta2), theta2, theta2)
% 
% dr6_eqn = diff(r6_eqn); % dr6 represented as diff(r6(theta2), theta2)
%% Calculate Values %%

theta3_array = ones(20,20);

for s = 0:20 % calculates values over the course of 20 seconds
%     theta2_num = t*dtheta2; % establish numerical value of theta2
    theta3_array(s) = subs(theta3_eqn, t, 1);
%     theta3_num = vpa(180 + solve(subs(theta3_eqn, t, s), theta3))
  
    %theta3_array(deg) = theta3_num
end

theta3_array

% theta3 eqn substitution with a numerical theta 2. 180 is a correction factor for the atan func

% dtheta3_num = subs(dtheta3_eqn, theta2, theta2_num);
% 
% ddtheta3_num = subs(ddtheta3_eqn, theta2, theta2_num);
% 
% r4_num = subs(r4_eqn, theta2, theta2_num);
% 
% dr4_num = subs(dr4_eqn, theta2, theta2_num);


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

dtheta2 = -15; % theta2 dot
ddtheta2 = 0; % theta2 doble-dot - second derivative

rho = % ENTER YOUR CODE HERE %; % density, gr/cm3
d = % ENTER YOUR CODE HERE %; % diameter, cm

m2 = % ENTER YOUR CODE HERE % ; % link 2, o2a2 kg
I_G4 = % ENTER YOUR CODE HERE %;
% and so on


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