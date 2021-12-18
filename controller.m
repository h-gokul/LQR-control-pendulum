%% Start
clear all;
clc;
close all;
%% Control of a two inverted pendulums on a cart
%% define variables 
syms dX X F M m1 m2 l1 l2 g;
syms x dx ddx theta1 dtheta1 ddtheta1 theta2 dtheta2 ddtheta2;

%% define state variables 
X = [x dx theta1 dtheta1 theta2 dtheta2];
constants = [M m1 m2 l1 l2 g];
u = [F];

%% Non-Linear system differential equations.

% double derivative of x (cart acceleration)
ddx_numerator = F - m1*g*cos(theta1)*sin(theta1) - - m2*g*cos(theta2)*sin(theta2) - m1*l1*(dtheta1^2)*sin(theta1) - m2*l2*(dtheta2^2)*sin(theta2);
ddx_denominator = M + m1*(sin(theta1))^2 + m2*(sin(theta2))^2;
ddx = ddx_numerator/ddx_denominator;

% double derivative of theta1
ddtheta1 = (ddx*cos(theta1) - g*sin(theta1))/l1;

% double derivative of theta2
ddtheta2 = (ddx*cos(theta2) - g*sin(theta2))/l2;

%  LHS of state eqn: X_dot
dX = [dx ddx dtheta1 ddtheta1 dtheta2 ddtheta2];

disp('====================================================================')
disp('*************************** Question 1 *****************************')
disp('Q1 a) Equations of motion of the Non-Linear system: ')
disp(['cart acceleration, ' ...
    'ddx ='])
disp(ddx)
disp(['pendulum 1 acceleration, ' ...
    'ddtheta1 ='])
disp(ddtheta1)
disp(['pendulum 2 acceleration, ' ...
    'ddtheta2 ='])
disp(ddtheta2)
disp('====================================================================')
%% Linearized State-space model
% (* defined along the equilibrium points)
% Equilibrium points:
X_e = [0 0 0 0 0 0];

% system matrix A can also be seen as the jacobian between dX and the states X
A = jacobian(dX, X);
A_origin = jacobian(dX, X);
B_origin = jacobian(dX, u);

A = subs(A, X, X_e); % defining A at equilibrium point
% input matrix B is nothing but the jacobian between dX and the inputs F
B = jacobian(dX, u);
B = subs(B, X, X_e); % defining B at equilibrium point
C = eye(size(A));
D = 0;
disp('====================================================================')
disp('Q1 b)State Space representation of our Linearized system: ')
disp('A Matrix: ')
disp(A)

disp('B Matrix: ')
disp(B)

disp('C Matrix: ')
disp(C)

disp('D Matrix: ')
disp(D)
disp('====================================================================')

%% Checking conditions on controllability using determinants
disp('====================================================================')
disp('Q1 c)Checking conditions on controllability : ')
controllability_mat = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
controllability_det = det(controllability_mat);
disp(controllability_det )
disp('====================================================================')

%% Substituting the values of constants
disp('====================================================================')
disp(['Q1 d) Substitute the values in the matrices, obtain LQR controller, ' ...
    'simulate response and certify local/global stability: '])
disp(" Substituting M =1000 Kg, m1=m2=100 Kg, l1=20m, l2=10m, g=9.8m/s")
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.8;

A = double(subs(A));
B = double(subs(B));

disp('A matrix:')
disp(A);
disp('B matrix:')
disp(B);

%% Check contrallability with values
disp("Rank of the controllability matrix : ")
disp(rank(ctrb(A,B)));
disp("It is full rank matrix, so the system at equilibrium is controllable.")

%% Check stability
openloop_poles = eigs(A);
disp('The Open loop poles are: ')
disp(openloop_poles);


%% LQR Design
disp('Desgining Linear Quadratic Regulator... ')
% Tuning Q
% Qx='1 ,' Qtheta1= '1 ,'Qtheta2='1 , 'R='1
% Qx='10 ,' Qtheta1= '100 ,'Qtheta2='200 , 'R='1
% Qx='100 ,' Qtheta1= '10000 ,'Qtheta2='20000 , 'R='1

% Tuning R
% Qx=1000; Qtheta1=1000000; Qtheta2=10000000; R_=0.0001; %ref

% Qx='2000 ,' Qtheta1= '2000000 ,'Qtheta2='1500000 , R=1
%Qx=2000; Qtheta1=1500000; Qtheta2=2000000; R_=0.1;

Qx=2000; Qtheta1=1500000; Qtheta2=2000000; R_=0.0001;
Q = zeros(6);
Q(1,1) = Qx;
Q(2,2) = 10;
Q(3,3) = Qtheta1;
Q(4,4) = 10;
Q(5,5) = Qtheta2;
Q(6,6) = 10;
disp('Q matrix : ')
disp(Q)

R = R_;
disp('R matrix : ')
disp(R)

K = lqr(A,B,Q,R);
disp('LQR gain matrix : ')
disp(K);


%% Closed Loop State-space system. 
disp('Closed Loop System Matrix (Ac): ')
Ac = A-B*K;

%% Check closed loop stability
disp('Closed Loop poles:')
closedloop_poles = eigs(Ac);
disp(closedloop_poles);
disp("all poles have real negative components, " + ...
    "closed loop the system is stable by the Lyapunov's indirect method.")
%% Simulation: 

% assuming initially at 5 degrees for theta1
x0 = 0.5; theta1_0 = deg2rad(10); theta2_0 = deg2rad(10);
init_state = [x0, 0, theta1_0, 0, theta2_0, 0];
openloop_state = ss(A, B, C, D);
closedloop_state = ss(Ac, B, C, D);

[y_op,t_op,x_op] = initial(openloop_state, init_state);
[y,t,x] = initial(closedloop_state, init_state);

% Unit step response
hold on
text = ['Unit step response of LQR Controller Qx=', Q(1,1) ,' Qtheta1= ', Q(3,3) ,'Qtheta2=', Q(3,3), 'R=', R];
plot(t(:,1),y(:,1),t(:,1),y(:,3), t(:,1),y(:,5));
legend x theta1 theta2
title(text)
hold off
% Q settles in 40 s tolerable response, transients are high tho.

disp('====================================================================')
disp('************************ end of Question 1 *************************')

disp('====================================================================')
disp('*************************** Question 2 *****************************')
disp('Q2 e) Check observability conditions: ')

%% Observability Check
disp('Checking Observability condition .. ')
% Only cart motion: x
C1 = [1 0 0 0 0 0];
disp(['Rank with only x(t)', num2str(rank(obsv(A,C1)))]);  
% rank=6 - observable

% For only pendulum motion: theta1, theta2
C2 = [0 0 1 0 0 0; 
      0 0 0 0 1 0];
disp(['Rank with theta1, theta2: ',num2str(rank(obsv(A,C2)))]);  
% rank=4 - not observable

% For cart and one pendulum motion: x, theta2
C3 = [1 0 0 0 0 0; 
      0 0 0 0 1 0];
disp(['Rank with x, theta2: ',num2str(rank(obsv(A,C3)))]);  
% rank=6 - observable

% For cart and 2x pendulum motion: x, theta1, theta2
C4 = [1 0 0 0 0 0; 
     0 0 1 0 0 0; 
     0 0 0 0 1 0];
disp(['Rank with x, theta1, theta2: ',num2str(rank(obsv(A,C4)))]);  
% rank=6 - observable
disp('====================================================================')

disp(' Qn f) Designing Luenberger Observer.. ')
%% Luenberger Observer
% 1) Find the smallest eigen value/pole in the closed loop system.

% 2) Set the desired pole to farther left of the left half plane 
% than the smallest pole. 
% 3) This choice of desired poles will result in faster convergence to 
% stability.

disp("##### For the Linearized System :")
disp("The real part of leftmost pole in the LQR closed loop system is")
disp(real(max(closedloop_poles))) % function takes absolute max idk why

% minimum real part was -2.76, lets choose desired poles lesser than that
disp("Lets choose the desired observer poles to be lesser than that")
desired_observer_poles= [-3.85 -4.75, -3.5, -3.4, -3.15, -2.9];
disp(desired_observer_poles);


% Observer pole placement and system design
% Case 1:  
L1 = place(A', C1', desired_observer_poles)';
A_L1 = [Ac (B*K)
    zeros(6,6) (A-L1*C1)];
B_L1 = [B 
        zeros(size(B))];
C_L1 = [C1 zeros(size(C1))];
D_L1 = 0;
observed_stateC1 = ss(A_L1, B_L1, C_L1, D_L1);

% Case 3:  
L3 = place(A', C3', desired_observer_poles)';
A_L3 = [Ac (B*K)
    zeros(6,6) (A-L3*C3)];
B_L3 = [B 
        zeros(size(B))];
C_L3 = [C3 zeros(size(C3))];
D_L3 = 0;
observed_stateC3 = ss(A_L3, B_L3, C_L3, D_L3);
% Case 4:  
L4 = place(A', C4', desired_observer_poles)';
A_L4 = [Ac (B*K)
    zeros(6,6) (A-L4*C4)];
B_L4 = [B 
        zeros(size(B))];
C_L4 = [C4 zeros(size(C4))];
D_L4 = 0;
observed_stateC4 = ss(A_L4, B_L4, C_L4, D_L4);

% Visualize Step input response: 
figure(1)
step(observed_stateC1)
title('Observability for Case1: Only cart motion: x(t)')
xlabel('Time')
ylabel('State')
% Visualize Step input response: 
figure(2)
step(observed_stateC3)
title('Observability for Case 3: cart and pendulum2 motion: x(t), theta2')
xlabel('Time')
ylabel('State')
% Visualize Step input response: 
figure(3)
step(observed_stateC4)
title('Observability for Case 4: cart and pendulum2 motion: x(t), theta2, theta3')
xlabel('Time')
ylabel('State')


%% LQG Output state feedback controller design
%% Optimal state estimator design:
%  we chose C1, the smallest observable output vector 
% [1 0 0 0 0 0]; 
D_out = zeros(size(C,1), size(B,2));

% Process and measurement noise 
w= 0.001*eye(size(A));  
v = 0.0001;
% optimal Kalman filter gain
Kf = lqe(A, w, C1, w, v);  
% observed state space system
A_kf = A - Kf*C1;
B_kf = [B Kf];
C_kf = eye(size(A));
D_kf = zeros(size(B_kf));
 
sysKF = ss(A_kf, B_kf, C_kf, D_kf);
 
%% Linear Quadratic Gaussion Controller for the Linearized System.
% LQG System = LQR + Kalman Filter
disp('====================================================================')
disp(' Qn f) Designing LQG Controller.. ')

A_lqg = [Ac, B*K;
         zeros(size(A)), (A*Kf*C1)];
B_lqg = [B; zeros(size(B))];

C_lqg = [C1 zeros(size(C1))];

closedloopsys_lqg = ss(A_lqg, B_lqg, C_lqg, 0);
X0 = [init_state init_state];
dt = 0.01; t_sim = dt:dt:100; u_in = ones(size(t_sim));
[y_lqg, t, x_lqg] = lsim(closedloopsys_lqg , u_in, t_sim, X0 );

% hold on
% figure(4)
% plot(t_sim, y_lqg(:, 1))
% xlabel('Cart position')
% ylabel('Time sec')
% hold off

% plot both response and estimated state trajectories of all components
% x_lqg contains the states and the error trajectories.
n_states = 6;
errors = x_lqg(:, n_states +1:end);
output_states = x_lqg(:, 1:n_states );
%estimated_states= output_states - errors;

% responses
x_response = output_states (:,1);
theta1_response = output_states(:,3);
theta2_response= output_states(:,5);
% estimates
% x_estimates= estimated_states(:,1);
% theta1_estimates= estimated_states(:,3);
% theta2_estimates= estimated_states(:,5);


hold on
subplot(3,1,1);
plot(t,x_response ,'-r')
legend x-response 
title 'LQG Controller: cart motion ';
% 
subplot(3,1,2);
plot(t,theta1_response,'-r')
legend theta1-response 
title 'LQG Controller: pendulum 1 motion';
% 
subplot(3,1,3);
plot(t,theta2_response,'-r')
legend theta2-response 
title 'LQG Controller: pendulum 2 motion';
disp("Eigen values of final LQR system")
disp(eig(A_lqg))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp("##### [TODO Qn f] Original Non linear system:")
% A_origin= subs(A_origin, constants, [M, m1, m2, l1, l2, g]) ;
% B_origin = subs(B_origin, constants, [M, m1, m2, l1, l2, g]);
% disp(A_origin)
% K_origin = lqr(A_origin,B_origin,Q,R); % choosing the same Q and R
% Ac_origin = A_origin-B_origin*K_origin;
% % Case 1:  
% L1_origin = place(A_origin', C1', desired_observer_poles)';
% A_L1_origin = [Ac_origin (B_origin*K)
%     zeros(6,6) (A_origin-L1_origin*C1)];
% B_L1_origin = [B_origin 
%         zeros(size(B_origin))];
% C_L1_origin = [C1 zeros(size(C1))];
% D_L1 = 0;
% observed_stateC1_origin = ss(A_L1_origin, B_L1_origin, C_L1_origin, D_L1);
% figure(1)
% step(observed_stateC1_origin)
% title('Observability for Case1: Only cart motion: x(t)')
% xlabel('Time')
% ylabel('State')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [UNCOMMENT BELOW]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Choose Cx with x =1,3,4 for Case1,Case3 or Case4 
% % estimate states using the luenberger observer for chosen case
% openloop_stateCx = ss(A, B, Cx, D);
% Lx = place(A', Cx', desired_observer_poles)';
% A_obCx = A-Lx*Cx; B_obCx = [B Lx]; 
% C_obCx = eye(size(A));D_obCx = zeros(size(B_obCx));
% 
% observed_stateCx = ss(A_obCx, B_obCx, C_obCx, D_obCx);
% dt = 0.01; t_sim = dt:dt:50; u_in = ones(size(t_sim));
% [y_openloop, t_] = lsim(openloop_state, u_in, t_sim);
% [y_openloopCx, t_] = lsim(openloop_stateCx, u_in, t_sim);
% [X_estimated, t_] = lsim(observed_stateCx, [u_in; y_openloopCx'], t_sim);
% % grid on
% % plot(t_,y_openloop) % unstable
% % plot(t_,y_openloopCx) % unstable
% 
% % hold on
% % subplot(3,1,1);
% % title 'Luenberger Observer: Tracking - x'
% % plot(t_(:,1), y_openloop(:,1),'k', t_(:,1), X_estimated(:,1), 'g--', 'Linewidth',1);
% % legend x estimated-x
% % hold off
% % 
% % hold on
% % subplot(3,1,2)
% % title 'Luenberger Observer: Tracking - theta 1'
% % plot(t_(:,1), y_openloop(:,3),'k', t_(:,1), X_estimated(:,3), 'g--', 'Linewidth',1);
% % legend theta1 estimate theta1
% % hold off
% % 
% % hold on
% % subplot(3,1,3)
% % title 'Luenberger Observer: Tracking - theta 2'
% % plot(t_(:,1), y_openloop(:,5),'k', t_(:,1), X_estimated(:,5), 'g--', 'Linewidth',1);
% % legend theta1 estimate theta1
% % hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% 
% disp('****End****')
% 
