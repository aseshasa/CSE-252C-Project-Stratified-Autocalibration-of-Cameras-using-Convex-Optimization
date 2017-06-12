function data = branchandbound(func, costfunc, data, tolerance)
f = inf;
for i = 1:40
% 	disp(i);
	[data, newb] = branch(data);
	try
		data = bound(func, data, newb);
	catch ME
		if length(ME.message) == 95 || length(ME.message) == 88
			break;
		else
			rethrow(ME);
		end
	end
	data = refine(costfunc, data);
	if abs(min(data.phi)-min(data.phi_lb)) <= tolerance
		break;
	elseif abs((min(data.phi)-f)/min(data.phi)) <= tolerance
		break;
	else
		f = min(data.phi);
	end
end
end

function data = bound(func, data, newb)
for i = newb
	[data.phi(i), data.phi_lb(i), data.qstar(:,i)] = func(data.l(:,i), data.u(:,i));
end
end

function [data, newb] = branch(data)
[~,idx1] = min(data.phi);
[~,idx2] = min(data.phi_lb);
idx = sort(unique([idx1,idx2]));
l1 = zeros(size(data.l,1),length(idx));
u1 = l1;
v1 = l1;
p1 = zeros(1,length(idx));
plb1 = p1;
for i = 1:length(idx)
	[~,dim] = max(abs(data.u(:,idx(i))-data.l(:,idx(i))));
	l1(:,2*i-1:2*i) = [data.l(:,idx(i)), data.l(:,idx(i))];
	l1(dim,2*i) = (data.l(dim,idx(i))+data.u(dim,idx(i)))/2;
	u1(:,2*i-1:2*i) = [data.u(:,idx(i)), data.u(:,idx(i))];
	u1(dim,2*i-1) = (data.l(dim,idx(i))+data.u(dim,idx(i)))/2;
	v1(:,2*i-1:2*i) = [data.qstar(:,idx(i)), data.qstar(:,idx(i))];
	p1(:,2*i-1:2*i) = [data.phi(:,idx(i)), data.phi(:,idx(i))];
	plb1(:,2*i-1:2*i) = [data.phi_lb(:,idx(i)), data.phi_lb(:,idx(i))];
	if v1(dim,2*i) >= l1(dim,2*i)
		p1(:,2*i-1:2*i) = [inf, data.phi(:,idx(i))];
		plb1(:,2*i-1:2*i) = [inf, data.phi_lb(:,idx(i))];
		isref(:,2*i-1:2*i) = [false, data.isref(:,idx(i))];
	else
		p1(:,2*i-1:2*i) = [data.phi(:,idx(i)), inf];
		plb1(:,2*i-1:2*i) = [data.phi_lb(:,idx(i)), inf];
		isref(:,2*i-1:2*i) = [false, data.isref(:,idx(i))];
	end
end
data.l(:,idx) = []; data.l = [data.l, l1];
data.u(:,idx) = []; data.u = [data.u, u1];
data.qstar(:,idx) = []; data.qstar = [data.qstar, v1];
data.phi(idx) = []; data.phi = [data.phi, p1];
data.phi_lb(idx) = []; data.phi_lb = [data.phi_lb, plb1];
data.isref(idx) = []; data.isref = [data.isref, isref];
newb = size(data.l,2)-2*length(idx)+1:size(data.l,2);
end