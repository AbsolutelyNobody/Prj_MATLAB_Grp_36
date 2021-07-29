%% initial parameters

%link lengths, units of metres
r1 = 7.8*10^(-2); % o2o3
r2 = 2.5*10^(-2); % o2a2
r3 = 13.8*10^(-2); % o3B
r5 = 4.75*10^(-2); % BC
r7 = 17.1*10^(-2); % o4o3

syms theta2 theta3 r4 theta5 r6;

% link 2 movement
dtheta2 = 1800; % deg/s from 300 rpm
ddtheta2 = 0; % constant angular velocity

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, using LCEs
% Found using link closure equations as detailed in the report
theta3_eqn = theta3 == atand((-sin(theta2))/(3.12 + cosd(theta2)));
r4_eqn = r4 == (2.5*sind(theta2)) / (sind(theta3));

theta5_eqn = theta5 == acosd((r7 - r3*cosd(180-theta3)) / (r5)) + 180;
r6_eqn = r6 == r3*sind(180-theta3) - r5*sin(theta5-180);

% Hint: Check if the angle needs to be adjusted to its true value
% Hint: Check this for all other angles too

%% Derivative equations of kinematic vars (d/dt) 
syms dtheta3 ddtheta3 dr4 dtheta5 ddtheta5 dr6

dtheta3_eqn = dtheta3 == diff(theta3_eqn);
ddtheta3_eqn = ddtheta3 == diff(dtheta3_eqn);

dr4_eqn = dr4 == diff(r4_eqn);

dtheta5_eqn = dtheta5 == diff(theta5_eqn);
ddtheta5_eqn = ddtheta5 == diff(dtheta5_eqn);

dr6_eqn = dr6 == diff(r6_eqn);
%% Calculate Values %%
theta2_num = 89;
vpa(subs(theta3_eqn, theta2, theta2_num))

%% Plot vars;

% Plot all desired deliverables. 
% 
% figure (1)
% plot(theta2,theta3_eqn)
% grid on;
% title('$\theta_3$ vs $\theta_2$', 'Interpreter','latex')
% xlabel('\theta_2   unit: degree')
% ylabel('\theta_3   unit: degree')
%  