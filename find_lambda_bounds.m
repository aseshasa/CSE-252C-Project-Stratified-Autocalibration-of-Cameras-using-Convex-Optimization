function [ L_bound, U_bound ] = find_lambda_bounds( H, l_mat, u_mat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_views = size(H,3);
h3 = H(:,3,:);
omega_0 = (l_mat + u_mat)/2;
L_bound = zeros(n_views,1);
U_bound = zeros(n_views,1);
opts = optimoptions(@fmincon, 'Display', 'off');
for i = 1:n_views
    h3_loc = h3(:,:,i);
    fun_min = @(x)h3_loc'*x*h3_loc;
    fun_max = @(x)-h3_loc'*x*h3_loc;
    l_b = max(l_mat, 0);
    A = [];
    b = [];
    Aeq = zeros(9);
	Aeq(9,9) = 1.0;
    beq = [0,0,0,0,0,0,0,0,1.0];
    %TODO: Normalization
    L_bound(i) = fun_min(fmincon(fun_max,omega_0, A,b,Aeq,beq,l_b, u_mat,{}, opts));
    U_bound(i) = fun_min(fmincon(fun_min,omega_0, A,b,Aeq,beq,l_b, u_mat,{}, opts));
end
L_bound = 1./L_bound;
U_bound = 1./U_bound;
end

