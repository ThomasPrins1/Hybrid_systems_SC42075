%% Variables
x_eff = 11.9; % mm
xr_eff = 16;  % mm
x_max = 40;   % mm
N_grind_max = 10;
psi = 0.98;
psi_real = 1.09;
gamma0 = 0;
gamma1 = 1;
gamma2 = 3;

% max variables:
j_max = 5; % number of sections
q_max = 3; % number of splits in a section
iter_max = 100; % number of iterations
%track_max = 70;
N_p = 10; % prediction horizon (in months?)

grp_nr = 4;
x_con = zeros(j_max,iter_max);
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
delta = np.zeros(2,1);

% cost_function
%J = J_deg + lapda*J_maint;
% J_deg = sum_j(sum_l(x_j_con)) summation of the length of degradation spots
% J_maint = sum_j(sum_l(sum_q(gamma*I))) 
% where I is an indication function to check if k+l-1 = a_q, its 1 if true,
% 0 if false
% gamma is given in the assignment
% l is the prediction window bounded iteration
% dynamics:

%% Functions
function out = predictTrackDeg(track,weightsA,weightsB)
    if track<30
        out = weightsA(1)*track+weightsB(1);
    elseif track<50
        out = weightsA(2)*track+weightsB(2);
    elseif track<70
        out = weightsA(3)*track+weightsB(3);
    else
        out = 0;
    end
end

function out = predictGrindingEffect(track,psi_f,x_eff_f)
    if track <= x_eff_f
        out = 0; % no grinding possible
    else
        out = psi_f*(track-x_eff_f); % grinding possible
    end
end

%% ! this part is no longer needed I believe!
for j = 1:j_max
    if x_con(j,1)<30
        %x_con_split(j,1) = x_con(j,1);
        f_degrade(j) = A(1)*x_con(j,1)+B(1);
    elseif x_con(j,1)<50
        %x_con_split(j,2) = x_con(j,1);
        f_degrade(j) = A(2)*x_con(j,1)+B(2);
    elseif x_con(j,1)<70
        %x_con_split(j,3) = x_con(j,1);
        f_degrade(j) = A(3)*x_con(j,1)+B(3);
    end
    % maybe we need this instead?
    %for q = 1:q_max
        %f_degrade(j,q) = A(q)*x_con_split(j,q)+B(q);
    %end
    if x_con(j,1) <= x_eff
        f_grind = 0; % no grinding possible
    else
        f_grind = psi*(x_con(j,1)-x_eff); % grinding possible
    end
% end of initialization
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
            f_degrade(j) = predictTrackDeg(x_con(j,k),A,B)
            x_con(j,k+1) = f_degrade(j);
        elseif u(j) == 1
            f_grind(j) = predictGrindingEffect(x_con(j,k),psi,x_eff)
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
    for l = 1:N_p
        % here we make our prediction based on l future iterations
        
    end
end
% l is the prediction horizon to look into the future
% our system is now a PWA with:
% x(k+1) = A_i*x(k) + B_i*u(k) + f_i
% y(k+1) = C_i*x(k)

% MLD part:
% first we check our states:
% x_con can be in 3 seperate states, so this requires 2 auxillary variables delta(1) and (2)
% [0 0] = state 0
% [1 0] = state 1
% [1 1] = state 2
% [0 1] = no state, make sure to exclude (by putting to 0?)
% also done here: https://www.researchgate.net/publication/262220797_MLD_Systems_Modeling_and_Control

%% Plotting
x_con_split;
f_degrade;