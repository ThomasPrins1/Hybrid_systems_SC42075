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

grp_nr = 4;

x_con = [20+0.5*(abs(grp_nr-15));
     22+0.5*(abs(grp_nr-15));
     25+0.5*(abs(grp_nr-15));
     27+0.5*(abs(grp_nr-15));
     17+0.5*(abs(grp_nr-15))];
x_aux = [7;8;7;7;8];
x0 = x_con;
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
j_max = 5; % number of sections
q_max = 3; % number of splits in a section
iter_max = 100; % number of iterations
u = zeros(1,j_max);
x_con_split = zeros(q_max,j_max);
f_degrade = zeros(q_max,j_max);

% cost_function
J = J_deg + lapda*J_maint;
% J_deg = sum_j(sum_l(x_j_con)) summation of the length of degradation spots
% J_maint = sum_j(sum_l(sum_q(gamma*I))) 
% where I is an indication function to check if k+l-1 = a_q, its 1 if true,
% 0 if false
% gamma is given in the assignment
% l is the prediction window bounded iteration
% dynamics:

%% Functions

%% Main code
% First part splits the x_con per section into degradation classes:
% 1->(0,30)
% 2->(30,50)
% 3->(50,70)
% f_degrade predicts the amount of degradation (squat size) in a section 
% after a month based on given weights in vector A, B and the split x_con
for j=1:j_max
    if x_con(j)<30
        x_con_split(1,j) = x_con(j);
    elseif x_con(j)<50
        x_con_split(2,j) = x_con(j);
    elseif (x_con(j)<70)
        x_con_split(3,j) = x_con(j);
    end
    for q=1:q_max
        f_degrade(q,j) = A(q)*x_con_split(q,j)+B(q);
    end
    if x_j_con <= x_eff
        f_grind = 0; % no grinding possible
    else
        f_grind = psi*(x_con-x_eff); % grinding possible
    end

end
for j=1:j_max
    for k = 1:iter_max % k is iterations so need to be reworked into other var
        % set degradation dynamics for each different value of u
        if u(j,k)== 0
            x_con(j,k) = f_degrade;
        elseif u(j,k) == 1
            x_con(j,k) = f_grind;
        elseif u(j,k) == 2 % replace
            x_con(j,k) = 0;
        end
        % set auxilary dynamics
        if u(j,k)== 0
            x_aux(j,k+1) = x_aux(j,k);
        elseif u(j,k) == 1
            x_aux(j,k+1) = x_aux(j,k)+1;
        elseif u(j,k) == 2 % replace
            x_aux(j,k+1) = 0;
        end
    end
end
x_con_split
f_degrade
%% Plotting