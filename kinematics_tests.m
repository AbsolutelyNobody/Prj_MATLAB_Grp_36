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
theta2(t) = t*dtheta2; % theta2 value 

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, using LCEs
% Found using link closure equations as detailed in the report

theta3(t) = atand((r2*sind(theta2(t)))/(r2*cosd(theta2(t))-r1))+180

r4(t)=(r2*cosd(theta2(t))-r1)/cosd(theta3(t));

theta5(t) = acosd((r7+r3*cosd(theta3(t))) / r5);

r6(t) = r5*sind(theta5(t)) - r3*sind(theta3(t));

% Hint: Check if the angle needs to be adjusted to its true value
% Hint: Check this for all other angles too

%% Derivative equations of kinematic vars (d/dt) 
syms dtheta3(t) ddtheta3(t) dr4(t) dtheta5(t) ddtheta5(t) dr6(t) ddr6(t)

dtheta3(t) = diff(theta3(t), t); % differentiates theta3_eqn with respect to t
ddtheta3(t) = diff(dtheta3(t), t); % differentiates dtheta3_eqn with respect to t, similar processes for the rest of the variables

dr4(t) = diff(r4(t), t);
ddr4(t) = diff(dr4(t), t);

dtheta5(t) = diff(theta5(t), t);
ddtheta5(t) = diff(dtheta5(t), t);

dr6(t) = diff(r6(t), t);
ddr6(t) = diff(dr6, t);

%% Calculate Values %%

theta3_array = ones(20, 1);

for s = 1:20 % calculates values over the course of 20 seconds

%     theta2_num = t*dtheta2; % establish numerical value of theta2
    theta3_array(s) = subs(theta3(t), t, s);
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
