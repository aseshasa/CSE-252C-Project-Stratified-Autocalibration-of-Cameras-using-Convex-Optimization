clear all;
close all;
cvx_solver sedumi;
w = warning ('off','all');
p0 = zeros(3,1);
f10 = 1;
f20 = 1;
s0 = 0;
u0 = 0;
v0 = 0;
K0 = [f10, s0, u0; 0, f20, v0; 0, 0, 1];
H0 = [K0, zeros(3,1); p0'*K0, 1];
trials = 50;
deltap = zeros(1,trials);
deltaf = deltap;
deltauv = deltap;
deltas = deltap;
for i = 1:trials
	disp(i);
	anim = generateToyAnimation(0, 'isProj', true, 'nPoint', 100, 'nFrame', 20);
	sigma = std(std(std(anim.W)));
	anim.W = anim.W + 0*sigma*randn(size(anim.W));
	anim.P = multiTimes(anim.P, inv(H0), 1);
	anim.S = H0*[anim.S; ones(1,anim.nPoint)];
	anim.S = anim.S(1:3,:) ./ (ones(3,1)*anim.S(4,:));
	[HEye, Hqa, p] = affineupgrade(anim, 1e-5);
	[H, K] = metric_upgrade(anim, p, 1e-5);
	deltap(i) = norm(p-p0);
	deltaf(i) = abs(K(1,1)/f10-1) + abs(K(2,2)/f20-1);
	deltauv(i) = abs(K(1,3)-u0) + abs(K(2,3)-v0);
	deltas(i) = abs(K(1,2)-s0);
	if deltap(i) > 1e-3
		break;
	end
end
format longg;
disp('Average deltap = '); disp(mean(deltap));
disp('Average deltaf = '); disp(mean(deltaf));
disp('Average deltauv = '); disp(mean(deltauv));
disp('Average deltas = '); disp(mean(deltas));