function [H, K] = metric_upgrade(anim, pInf, tolerance)
data.l = [0.9, -0.1, -0.1,0.9,-0.1]';
data.u = [1.1, 0.1, 0.1,1.1,0.1]';

num_views = anim.nFrame;
H_inf = zeros(3,3,num_views);

for i = 1:num_views
	H_inf(:,:,i) = anim.P(:,1:3,i) - anim.P(:,4,i)*pInf';
end

hand1 = @(x)costfunc(x,H_inf);
hand2 = @(x,y)eqn28(H_inf, x, y);
data.qstar = (data.l+data.u)/2;
data.phi = -Inf;
data.phi_lb = -Inf;
data.isref = false;
data = refine(hand1, data);
data = branchandbound(hand2, hand1, data, tolerance);
% [data.phi, data.phi_lb, data.k] = eqn28(H_inf, data.l, data.u);
[~, idx] = min(data.phi);
k = data.qstar(:,idx);
K = [k(1:3)'; 0 k(4:5)'; 0 0 1];
H = [[K;-pInf'*K],[0;0;0;1]];
end

function [phi, phi_lb, k] = eqn28(H_inf, l, u)
[l_mat, u_mat] = prop_interval(l, u);
[L_lambda, U_lambda] = find_lambda_bounds(H_inf, l_mat,u_mat);
num_views = size(H_inf,3);
cvx_begin sdp
	variable omega(3,3) nonnegative semidefinite;
	variable lambda(num_views);
	variable v(num_views, 3,3);
	expression t(num_views);
	for i = 1:num_views
		t(i) = norm(omega - squeeze(H_inf(:,:,i)) * squeeze(v(i,:,:)) * squeeze(H_inf(:,:,i))', 'fro');
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
phi_lb = cvx_optval;
k = omega2k(omega);
phi = costfunc(k, H_inf);
end

function cost = costfunc(k,H_inf)
omega = k2omega(k);
nFrame=size(H_inf,3); nSample=size(k,2);
cost=zeros(1,nSample);
for i=1:nFrame
	HOmegaHt=multiTimes(H_inf(:,:,i),multiTimes(omega,H_inf(:,:,i)',1),1.2);
	cost=cost+sum(reshape(omega-bsxfun(@rdivide,HOmegaHt,HOmegaHt(3,3,:)),9,nSample).^2,1);
end
end

function k = omega2k(omega)
[U,S,V]=svd(omega); S(S<0)=0; omega2=U*S*V';
L=zeros(3,3);
omega2=omega2(end:-1:1,end:-1:1);
for j=1:3
	for i=j:3
		if i==j
			L(i,i)=sqrt(max(omega2(i,i)-sum(L(i,1:i-1).^2),0));
		else
			L(i,j)=(omega2(i,j)-L(i,1:j-1)*L(j,1:j-1)')/L(j,j);
		end
	end
end
L=L';
K=L(end:-1:1,end:-1:1)';
K=K/K(3,3);
k = [K(1,1:3)'; K(2,2:3)'];
end

function omega = k2omega(k)
K = zeros(3,3,size(k,2));
K(1,1,:)=k(1,:);
K(1,2,:)=k(2,:);
K(1,3,:)=k(3,:);
K(2,2,:)=k(4,:);
K(2,3,:)=k(5,:);
K(3,3,:)=1;
omega = multiTimes(K,K,2.2);
end