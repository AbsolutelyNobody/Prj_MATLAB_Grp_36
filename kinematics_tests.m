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
subs(theta3_eqn, t, 0.2)
%for s = 0:20 % calculates values over the course of 20 seconds
%     theta2_num = t*dtheta2; % establish numerical value of theta2
%     
%     theta3_num = vpa(180 + solve(subs(theta3_eqn, t, s), theta3))
  
    %theta3_array(deg) = theta3_num
%end

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
