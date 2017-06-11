function [H, K] = metric_upgrade(anim, pInf)
    l_k = [[0.9, -0.1, -0.1];[0,0.9,-0.1];[0,0,1]];
    u_k = [[1.1, 0.1, 0.1];[0,1.1,0.1];[0,0,1]];
    [l_mat, u_mat] = prop_interval(l_k, u_k);
    
    num_views = anim.nFrame;
    H_inf = zeros(num_views, 3,3);
    
    for i = 1:num_views
        H_inf(i,:,:) = anim.P(:,1:3,i) - anim.P(:,4,i)*pInf';
    end
    [L_lambda, U_lambda] = find_lambda_bounds(H_inf, l_mat,u_mat);

    cvx_begin sdp quiet
        variable omega(3,3) nonnegative semidefinite;
        variable lambda(num_views);
        variable v(num_views, 3,3);
        expression t(num_views);
        for i = 1:num_views
           t(i) = norm(omega - squeeze(H_inf(i,:,:)) * squeeze(v(i,:,:)) * squeeze(H_inf(i,:,:))', 'fro');
        end
        minimize(sum(pow_pos(t,2)));
        subject to
            l_mat(:) <= omega(:) <= u_mat(:);
            L_lambda(:) <= lambda(:) <= U_lambda(:);
            omega >= 0;
            omega(3,3) == 1;
            for i = 1:num_views
                for j = 1:3
                    for k = 1:3
                        v(i,j,k) <= U_lambda(i) * omega(j,k) + l_mat(j,k) * lambda(i) - U_lambda(i) * l_mat(j,k);
                        v(i,j,k) <= L_lambda(i) * omega(j,k) + u_mat(j,k) * lambda(i) - L_lambda(i) * u_mat(j,k);
                        v(i,j,k) >= L_lambda(i) * omega(j,k) + l_mat(j,k) * lambda(i) - L_lambda(i) * l_mat(j,k);
                        v(i,j,k) >= U_lambda(i) * omega(j,k) + u_mat(j,k) * lambda(i) - U_lambda(i) * u_mat(j,k);
                    end
                end
            end
    cvx_end
    
    
    K = chol(omega);
    H = [[K;-pInf'*K],[0;0;0;1]];
end