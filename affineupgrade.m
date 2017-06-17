function [HEye, Hqa, p] = affineupgrade(anim, tolerance)
[HEye, Hqa, l, u, P] = chiralineq(anim);

P11=squeeze(P(1,1,:));
P12=squeeze(P(1,2,:));
P13=squeeze(P(1,3,:));
P14=squeeze(P(1,4,:));
P21=squeeze(P(2,1,:));
P22=squeeze(P(2,2,:));
P23=squeeze(P(2,3,:));
P24=squeeze(P(2,4,:));
P31=squeeze(P(3,1,:));
P32=squeeze(P(3,2,:));
P33=squeeze(P(3,3,:));
P34=squeeze(P(3,4,:));
alpha = -[ P14, P24, P34, - P33 - P22 - P11 ];
beta = [ P12.*P24 + P34.*P13 - P33.*P14  - P14.*P22, P34.*P23 - P11.*P24 + ...
  P14.*P21 - P33.*P24, - P11.*P34 - P34.*P22 + P32.*P24 + P14.*P31, P11.*P33 + ...
  P11.*P22 - P12.*P21 + P33.*P22 - P32.*P23 - P13.*P31 ];
gamma = -[ P34.*P12.*P23 - P34.*P13.*P22 - P33.*P12.*P24 - P32.*P14.*P23 + ...
  P33.*P14.*P22 + P32.*P13.*P24, -P11.*P34.*P23 + P34.*P13.*P21 + P11.*P33.*P24 - ...
  P33.*P14.*P21 - P13.*P24.*P31 + P14.*P23.*P31, P11.*P34.*P22 - P34.*P12.*P21 - ...
  P11.*P32.*P24 + P32.*P14.*P21 + P12.*P24.*P31 - P14.*P22.*P31, P33.*P12.*P21 - ...
  P11.*P33.*P22 + P11.*P32.*P23 - P32.*P13.*P21 - P12.*P23.*P31 + P13.*P22.*P31 ];

aH = alpha*Hqa';
bH = beta*Hqa';
cH = gamma*Hqa';
dH = ones(size(aH,1),1)*Hqa(:,4)';
data.l = l(1:3); data.u = u(1:3);
data.qstar = (l(1:3)+u(1:3))/2;
data.phi = -Inf;
data.phi_lb = -Inf;
data.isref = false;
hand1 = @(x)costfunc(x,aH,bH,cH,dH);
hand2 = @(x,y)eqn21(aH, bH, cH, dH, x, y);
data = branchandbound(hand2, hand1, data, tolerance);
% [data.phi, data.phi_lb, data.qstar] = eqn21(aH, bH, cH, dH, data.l, data.u);
[~, idx] = min(data.phi);
p=Hqa'*[data.qstar(:,idx);1];
p=p(1:3)/p(4);
end

function [phi, phi_lb, v] = eqn21(aH, bH, cH, dH, l, u)
vbd = reshape([l, u; 1 1],1,4,2);
abd = sum(sort(bsxfun(@times, aH, vbd),3),2);
bbd = sum(sort(bsxfun(@times, bH, vbd),3),2);
cbd = sum(sort(bsxfun(@times, cH, vbd),3),2);
dbd = sum(sort(bsxfun(@times, dH, vbd),3),2);
a_l = abd(:,:,1);
a_u = abd(:,:,2);
b_l = bbd(:,:,1);
b_u = bbd(:,:,2);
c_l = cbd(:,:,1);
c_u = cbd(:,:,2);
d_l = dbd(:,:,1);
d_u = dbd(:,:,2);
m = size(aH,1);
cvx_begin sdp quiet
	variables v(3) r e f(m,1) g(m,1) tf1(m,1) tg1(m,1) ypf2(m,4) ypg2(m,4);
	minimize(r);
	subject to
		a = aH*[v;1];
		b = bH*[v;1];
		c = cH*[v;1];
		d = dH*[v;1];
		{d(1),e} <In> relax_x83(d(1), d_l(1), d_u(1), e);
		l <= v <= u;
		{(f - g), r, e} <In> rotated_lorentz(m);
		{f,c,a,tf1,ypf2} <In> x13y_relax(f, c, c_l, c_u, a, a_l, a_u, tf1, ypf2);
		{g,d,b,tg1,ypg2} <In> x13y_relax(g, d, d_l, d_u, b, b_l, b_u, tg1, ypg2);
cvx_end
cvx_clear;
phi_lb = r;
phi = costfunc(v, aH, bH, cH, dH);
end

% relaxations for z = x^(8/3)
function ret = relax_x83(x, xl, xu, z)
xl13 = nthroot(xl,3);
xu13 = nthroot(xu,3);
cvx_begin
	subject to
		z <= xl13.^8+(x-xl).*(xu13.^8-xl13.^8)./(xu-xl);
cvx_end
ret = {x, z};
end

% relaxations for z = x^(1/3)*y
function ret = x13y_relax(z, x, xl, xu, y, yl, yu, t, yp)
idx_case11 = xl>0;
idx_case12 = xu<0;
idx_case2 = (xl<=0) & (xu>=0);
cvx_begin
	subject to
		if sum(idx_case11)
			{x(idx_case11),y(idx_case11),z(idx_case11),yp(idx_case11,1)} <In> case1_relax(x(idx_case11), ...
				xl(idx_case11),xu(idx_case11),y(idx_case11),yl(idx_case11),yu(idx_case11),z(idx_case11),yp(idx_case11,1));
			{x(idx_case11),-y(idx_case11),-z(idx_case11),yp(idx_case11,2)} <In> case1_relax(x(idx_case11), ...
				xl(idx_case11),xu(idx_case11),-y(idx_case11),-yu(idx_case11),-yl(idx_case11),-z(idx_case11),yp(idx_case11,2));
		end
		if sum(idx_case12)
			{-x(idx_case12),-y(idx_case12),z(idx_case12),yp(idx_case12,3)} <In> case1_relax(-x(idx_case12), ...
				-xu(idx_case12),-xl(idx_case12),-y(idx_case12),-yu(idx_case12),-yl(idx_case12),z(idx_case12),yp(idx_case12,3));
			{-x(idx_case12),y(idx_case12),-z(idx_case12),yp(idx_case12,4)} <In> case1_relax(-x(idx_case12), ...
				-xu(idx_case12),-xl(idx_case12),y(idx_case12),yl(idx_case12),yu(idx_case12),-z(idx_case12),yp(idx_case12,4));
		end
		if sum(idx_case2)
			{x(idx_case2),y(idx_case2),z(idx_case2),t(idx_case2)} <In> case2_relax(x(idx_case2),xl(idx_case2), ...
				xu(idx_case2),y(idx_case2),yl(idx_case2),yu(idx_case2),z(idx_case2),t(idx_case2));
		end
cvx_end
ret = {z, x, y, t, yp};
end

% convex relaxation for z = x^(1/3)*y, xl > 0 or xu < 0
function ret = case1_relax(x, xl, xu, y, yl, yu, z, yp)
xu13 = nthroot(xu, 3);
xl13 = nthroot(xl, 3);
cvx_begin
	subject to
		lambda = (x-xl)./(xu-xl);
		z >= xl13.*yp + xu13.*(y - yp);
		(1 - lambda).*yl <= yp <= (1 - lambda).*yu;
		lambda.*yl <= y - yp <= lambda.*yu;
cvx_end
ret = {x, y, z, yp};
end

% relaxations for z = x^(1/3)*y, xl <= 0 <= xu 
function ret = case2_relax(x, xl, xu, y, yl, yu, z, t)
idxconv1 = -xu/8 > xl;
idxconv2 = ~idxconv1;
idxconc1 = -xl/8 < xu;
idxconc2 = ~idxconc1;
xl13 = nthroot(xl, 3);
xu13 = nthroot(xu, 3);
xum23 = xu13.^(-2);
mxl23 = (-xl13).^2;
xlm23 = xl13.^(-2);
xu23 = xu13.^2;
cvx_begin
	subject to
		if sum(idxconv1)
			t(idxconv1) >= (xum23(idxconv1) - 2/3./xu(idxconv1).*xl13(idxconv1)).*x(idxconv1) + ...
				2/3.*xl13(idxconv1);
			t(idxconv1) >= (x(idxconv1)+2*xl(idxconv1))./3./mxl23(idxconv1);
		end
		if sum(idxconc1)
			t(idxconc1) <= (xlm23(idxconc1) - 2/3./xl(idxconc1).*xu13(idxconc1)).*x(idxconc1) + ...
				2/3.*xu13(idxconc1);
			t(idxconc1) <= (x(idxconc1)+2*xu(idxconc1))./3./xu23(idxconc1);
		end
		if sum(idxconv2)
			t(idxconv2) >= ((xu13(idxconv2)-xl13(idxconv2)).*x(idxconv2) + (xu(idxconv2).* ...
				xl13(idxconv2)-xl(idxconv2).*xu13(idxconv2)))./(xu(idxconv2)-xl(idxconv2));
		end
		if sum(idxconc2)
			t(idxconc2) <= ((xu13(idxconc2)-xl13(idxconc2)).*x(idxconc2) + (xu(idxconc2).* ...
				xl13(idxconc2)-xl(idxconc2).*xu13(idxconc2)))./(xu(idxconc2)-xl(idxconc2));
		end
		{t, y, z} <In> bilin_relax(t, xl13, xu13, y, yl, yu, z);
cvx_end
ret = {x, y, z, t};
end

% relaxations for z = x*y;
function ret = bilin_relax(x, xl, xu, y, yl, yu, z)
cvx_begin
	subject to
		z >= xl.*y + yl.*x - xl.*yl;
		z >= xu.*y + yu.*x - xu.*yu;
		z <= xu.*y + yl.*x - xu.*yl;
		z <= xl.*y + yu.*x - xl.*yu;
cvx_end
ret = {x, y, z};
end

function [cost,grad]=costfunc(v, aH, bH, cH, dH)
v(4,:)=1;
a = aH*v;
b = bH*v;
c13 = nthroot(cH*v, 3);
d13=nthroot(dH*v,3);
cost=sum((c13.*a-bsxfun(@times,d13,b)).^2,1)./d13(1).^8;
end