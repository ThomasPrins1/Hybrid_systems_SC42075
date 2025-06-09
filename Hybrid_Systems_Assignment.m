%% Variables
x_eff = 11.9; % mm
xr_eff = 16;  % mm
x_max = 40;   % mm
N_grind_max = 10;
psi = 0.98;
psi_real = 1.09;
gamma1 = 0;
gamma2 = 1;
gamma3 = 3;
gamma = [gamma1,gamma2,gamma3];

% max variables:
j_max = 5; % number of sections
q_max = 3; % number of splits in a section
iter_max = 100; % number of iterations
N_p = 2; % prediction horizon (in months?)

grp_nr = 4;
x_con = zeros(j_max,iter_max+N_p);
x_con(:,1) = [20+0.5*(abs(grp_nr-15));
     22+0.5*(abs(grp_nr-15));
     25+0.5*(abs(grp_nr-15));
     27+0.5*(abs(grp_nr-15));
     17+0.5*(abs(grp_nr-15))];
x_aux = zeros(j_max,iter_max);
x_aux(:,1) = [7;8;7;7;8];

A_real = [1.0017;
          1.015;
          1.008];
B_real = [1.6959;
          2.2165];

A = [1.0075;
     1.0095;
     1.01];
B = [1.65+0.001*abs(grp_nr-15);
     2.07+0.001*abs(grp_nr-15);
     2.47+0.001*abs(grp_nr-15)];

% matrix definitions
u = zeros(j_max,1);
x_con_split = zeros(j_max,q_max);
f_degrade = zeros(j_max,1);
f_grind = zeros(j_max,1);

% MLD
delta = zeros(3,1);

% cost_function
% J = J_deg + lapda*J_maint;
% J_deg = sum_j(sum_l(x_j_con)) summation of the length of degradation spots
% J_maint = sum_j(sum_l(sum_q(gamma*I))) 
% where I is an indication function to check if k+l-1 = a_q, its 1 if true,
% 0 if false
% gamma is given in the assignment
% l is the prediction window bounded iteration
% dynamics:

%% Functions
function out = convexMILP(x_con,lapda,gamma,k,n,N_p,A,B,x_eff,psi)
    % need to save x_con in convex variable:
    x_pred = optimexpr(n, N_p + 1);  
    x_pred(:,1) = x_con(:,k);
    x0.delta_u = [1; 0; 0];   % Initial guess for delta_u
    x0.delta = [1; 0; 0];     % Initial guess for delta
    x0.z = [10; 40; 60];      % Initial guess for z (within bounds)
    x0.delta_grind = [0];
    x0.z_grind = [50];
    delta_u = optimvar('delta_u', 3, 1, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
    delta = optimvar('delta', 3, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
    delta_grind = optimvar('delta_grind', 1, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
    z = optimvar('z', 3, 'LowerBound', 0);
    z_grind = optimvar('z_grind', 1, 'LowerBound', 0);

    x_min = [0,30,50];
    x_max = [30,50,70];
    ineq_constraints = [];
    eq_constraints = [];
    eq_constraints = [eq_constraints; sum(delta_u) == 1];
    ineq_constraints = [ineq_constraints;
        -delta(1) + delta(3) <= 0;
        -delta(2) + delta(3) <= 0;
        delta(1) + delta(2) - delta(3) <= 1;
    ];
    
    for q = 1:3
        ineq_constraints = [ineq_constraints;
            z(q) <= x_max(q) * delta(q);
            -z(q) <= -x_min(q) * delta(q); % multiplied by -1 to change to <=
        ];
    end
    % u is also a binary value that we want to use inside our convex 
    % optimization, so this also requires a auxilary variable to decide if
    % MPC will choose to do maintenance, replace or do nothing
    for j = 1:n
        for l = 1:N_p
            track_prediction = predictTrackDeg(x_pred(j,k+l-1),A,B,delta,z,"MLD");
            [grinding_prediction,x_grind_max,x_grind_min] = predictGrindingEffect(x_pred(j,k+l-1),psi,x_eff,delta_grind,z_grind,"MLD");
            x_pred(j,k+l) = delta_u(1)*track_prediction + delta_u(2)*grinding_prediction + delta_u(3)*0;
            for q = 1:3
                ineq_constraints = [ineq_constraints;
                z(q) <= x_pred(j,k+l-1) - x_min(q) * (1 - delta(q));
                -z(q) <= -x_pred(j,k+l-1) + x_max(q) * (1 - delta(q)); % multiplied by -1 to change to <=
                ];
            end
            ineq_constraints = [ineq_constraints;
            z_grind <= x_pred(j,k+l-1) - x_grind_min * (1 - delta_grind);
            -z_grind <= -x_pred(j,k+l-1) + x_grind_max * (1 - delta_grind);
            ];
        end
    end
    ineq_constraints = [ineq_constraints;
        z_grind <= x_grind_max * delta_grind;
        -z_grind <= -x_grind_min * delta_grind;
        ]
    J_deg_pred = sum(x_pred);
    J_maint = 0;
    for j = 1:n
        for l = 1:N_p
            for q = 1:3
                add_JMaint = gamma(q) * delta_u(q);
                J_maint = J_maint + add_JMaint;
                % assigment says gamma(j,q), but i cannot see how gamma is different between tracks?
                % This would onyl be a small change anyway, but gamma
                % values are only given for q=1,q=2 and q=3, without
                % differences between j values
            end
        end
    end
    J_maint_lapda = lapda*J_maint;

    J = J_deg_pred + J_maint_lapda;

    prob = optimproblem;
    prob.Objective = J;
    prob.Constraints.eq = eq_constraints;
    prob.Constraints.ineq = ineq_constraints;
    options = optimoptions('intlinprog', 'Display', 'iter'); % Show solver output
    %int_vars = [delta_u; delta];
    % Solve the problem
    x_pred(:,1)
    [sol, fval, exitflag, output] = solve(prob, x0, 'Options', options, 'Solver', 'intlinprog');
    [~, u] = max(sol.delta_u);
    % output u for MPC given:
    out = u - 1; % Because MATLAB indices start at 1, so subtract 1 to get u ∈ {0,1,2}
end

function out = predictTrackDeg(track,weightsA,weightsB,delta,z,controller_type)
    if controller_type == "PWA"
        if track<30
            out = weightsA(1)*track+weightsB(1);
        elseif track<50
            out = weightsA(2)*track+weightsB(2);
        elseif track<70
            out = weightsA(3)*track+weightsB(3);
        else
            out = 0;
        end
    elseif controller_type == "MLD"
        % fully written out the following needs to hold:
        %dynamic1 = ((1-delta(1))*(1-delta(2))*weightsA(1)*track+weightsB(1));
        %dynamic2 = (delta(1)*(1-delta(2))*weightsA(2)*track+weightsB(2));
        %dynamic3 = ((1-delta(1))*delta(2)*weightsA(3)*track+weightsB(3));
        % This results in terms with delta1*delta2, which we set equal to
        % delta3 under the following conditions:
        % -delta1+delta3 <= 0
        % -delta1+delta3 <= 0
        % delta1 + delta2 - delta3 <= 1
        %delta(3) = delta(1)*delta(2); % This is held by constraints!
        %dynamic1 = (1+delta3-delta(1)-delta(2)*weightsA(1)*track+weightsB(1));
        %dynamic2 = ((delta(2)-delta(3))*weightsA(2)*track+weightsB(2));
        %dynamic3 = ((delta(1)-delta(3))*weightsA(3)*track+weightsB(3));
        % Then the next step will be to take out the delta1*x and delta2*x
        % and replace this with auxillary variables z1 and z2, this will
        % cause the function to be represented by z, delta, x all in
        % liniear form (similarly included in convex!):
        %z1 = delta(1)*track;
        %z2 = delta(2)*track;
        %z3 = delta(3)*track;
        %z = [z1,z2,z3]
        dynamic1A = -z(1)*weightsA(1)-z(2)*weightsA(1)+z(3)*weightsA(1)+track*weightsA(1);
        dynamic1B = -delta(1)*weightsB(1)-delta(2)*weightsB(1)+delta(3)*weightsB(1)+weightsB(1);
        dynamic2 = z(2)*weightsA(2) - z(3)*weightsA(2) + delta(2)*weightsB(2) - delta(3)*weightsB(2);
        dynamic3 = z(1)*weightsA(3) - z(3)*weightsA(3) + delta(1)*weightsB(3) - delta(3)*weightsB(3);
        % Now z also needs to hold for some conditions:
        % z <= max(x)*delta
        % z >= min(x)*delta
        % z <= x - min(x)(1-delta)
        % z >= x - max(x)(1-delta)
        out = dynamic1A + dynamic1B + dynamic2 + dynamic3;
        % 4th dynamic with delta1 = 1 && delta 2 equal to 1 will just be
        % equal to 0
    end
end

function [out,max,min] = predictGrindingEffect(track,psi_f,x_eff_f,delta_grind,z_grind,controller_type)
    max = 70
    min = x_eff_f
    if controller_type == "PWA"
        if track <= x_eff_f
            out = 0; % no grinding possible
        else
            out = psi_f*(track-x_eff_f); % grinding possible
        end
    elseif controller_type == "MLD"
        
        %[out,max] = delta_grind*(psi_f*(track-x_eff_f));
        out = psi_f*(z_grind - delta_grind*x_eff_f);
    end
end

%% Main code
% First part splits the x_con per section into degradation classes:
% 1->(0,30)
% 2->(30,50)
% 3->(50,70)
% f_degrade predicts the amount of degradation (squat size) in a section 
% after a month based on given weights in vector A, B and the split x_con

for k = 1:iter_max % k is iterations so need to be reworked into other var
    for j = 1:j_max
        % set degradation dynamics for each different value of u
        if u(j)== 0 
            f_degrade(j) = predictTrackDeg(x_con(j,k),A,B,delta,[0,0,0],"PWA")
            x_con(j,k+1) = f_degrade(j);
        elseif u(j) == 1
            [f_grind(j),~,~] = predictGrindingEffect(x_con(j,k),psi,x_eff,delta,0,"PWA") % constant so nothing needed
            x_con(j,k+1) = f_grind(j);
        elseif u(j) == 2 % replace
            x_con(j,k+1) = 0;
        end
        % set auxilary dynamics, which count the amount of grindings done.
        if u(j)== 0 % no action
            x_aux(j,k+1) = x_aux(j,k);
        elseif u(j) == 1 && x_aux(j,k) < N_grind_max % grind
            x_aux(j,k+1) = x_aux(j,k)+1;
        elseif u(j) == 2 % replace
            x_aux(j,k+1) = 0;
        end
    end

    % MPC part:
    % lapda = linspace(600,800,50)
    % for r = 1:51 ...
    lapda = 600; % needs to be between 600-800
    u(k) = convexMILP(x_con,lapda,gamma,k,j_max,N_p,A,B,x_eff,psi)
    %predictTrackDeg(x_con(j,k),A,B,delta,"MLD")
end
% l is the prediction horizon to look into the future
% our system is now a PWA with:
% x(k+1) = A_i*x(k) + B_i*u(k) + f_i
% y(k+1) = C_i*x(k)

% MLD part:
% first we check our states:
% x_con can be in 3 seperate states, so this requires 2 auxillary variables delta(1) and (2)
% smaller then 30, smaller then 50, smaller then 70:
% [0 0] = state 0 -> 30 gives A1*x+B1
% [1 0] = state 1 -> 50 gives A2*x+B2
% [1 1] = state 2 -> 70 gives A3*x+B3
% [0 1] = no state, make sure to exclude (by putting to 0?)
% also done here: https://www.researchgate.net/publication/262220797_MLD_Systems_Modeling_and_Control
% We create an extended state space: x(k) = [x(k), delta1, delta2]
% if delta1 = 0 & delta2 = 0: A1*x + B1*u + 0 + 0
% if delta1 = 1 & delta2 = 0: A1 + A2 + B1*u + B2 + 0 = A2 + B2*u
% if delta1 = 0 & delta2 = 1: A1 + A3 + B1*u + B3 + 0 = A3 + B3*u
% A123 + B1 A2+B2 A3+B3

% Furthermore, our control input is also discrete.
% 
%% Plotting
x_con_split;
f_degrade;