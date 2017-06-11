function data = refine(func, data)
[~,idx1] = min(data.phi_lb);
idx2 = find(data.phi_lb < min(data.phi_lb(data.isref)) & ~data.isref);
if ~data.isref(idx1)
	idx = unique([idx2, idx1]);
else
	idx = idx2;
end
for i = idx
	temp = min([max([data.qstar(:,i),data.l(:,i)],[],2),data.u(:,i)],[],2);
	randval = rand(size(data.l,1),10000) .* ((data.u(:,i)-data.l(:,i))*ones(1,10000));
	randq = [temp, randval + data.l(:,i)*ones(1,10000)];
	randcost = func(randq);
	[~,ind] = min(randcost);
	temp = randq(:,ind);
	temp = fmincon(func,temp,[],[],[],[],data.l(:,i),data.u(:,i),[],optimset(...
    'GradObj','off','Hessian','off','Algorithm','active-set','Display','off'));
	temp = min([max([temp,data.l(:,i)],[],2),data.u(:,i)],[],2);
	tempcost = func(temp);
	if tempcost<data.phi(i)
		data.phi(i) = tempcost; data.qstar(:,i) = temp;
	end
	data.isref(i)=true;
end
bad = find(data.phi_lb > min(data.phi));

if ~isempty(bad) && length(bad)~=length(data.phi)
  data.l(:,bad)=[]; data.u(:,bad)=[]; data.qstar(:,bad)=[];
  data.phi_lb(bad)=[]; data.phi(bad)=[];
  data.isref(bad)=[];
end