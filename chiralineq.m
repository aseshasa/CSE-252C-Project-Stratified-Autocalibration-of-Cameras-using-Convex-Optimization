function [HEye, Hqa, l, u, P] = chiralineq(anim)
[anim,HEye]=anim.setFirstPRtToId();
x = [anim.W; ones(1, anim.nPoint, anim.nFrame)];
X = [anim.S; ones(1, anim.nPoint)];
P = anim.P;
nFrame = anim.nFrame;
nPoint = anim.nPoint;

% find wi's
w = zeros(1, nPoint, nFrame);
for i = 1:nFrame
	xhat = P(:,:,i)*X;
	w(:,:,i) = mean(xhat./x(:,:,i), 1);
end

%change signs of P's, X's
if sum(sign(w(:,:,1))) < 0
	X = -X;
	w = -w;
end
for i = 1:nFrame
	if sum(sign(w(:,:,i))) < 0
		P(:,:,i) = -P(:,:,i);
		w(:,:,i) = -w(:,:,i);
	end
end
for i = 1:nPoint
	if sum(sign(w(:,i,:))) < 0
		X(:,i,:) = -X(:,i,:);
		w(:,i,:) = -w(:,i,:);
	end
end

%find camera centres using nullspace of P
C = zeros(4, nFrame);
for i = 1:nFrame
	C(:,i) = null(P(:,:,i));
	C(:,i) = C(:,i)/C(4,i);
end

%alg 21.1: Find H
for delta = [-1, 1]
	cvx_begin quiet
		variables v(4) d;
		maximize(d);
		subject to
			X'*v >= d;
			delta*C'*v >= d;
			-1 <= v <= 1;
	cvx_end
	cvx_clear;
	if d >= 0 && any(abs(v)>1e-5)
		break;
	end
end
Hq = [[eye(3), zeros(3,1)]; v'];
Hq(1,1) = delta/Hq(4,4);

%alg 21.2: Find Hqa, v_lim
HqX = Hq*X;
HqX = HqX(1:3,:)./(ones(3,1)*HqX(4,:));
mu = mean(HqX, 2);
HqXc = HqX - repmat(mu, 1, nPoint);
[U, S] = eig(HqXc*HqXc');
S = diag(1./sqrt(diag(S)));
K = S*U';
if det(U) < 0
	K = -K;
end
Ha = [K -K*mu; 0 0 0 1];
Hqa = Ha*Hq;
Xqa = Hqa*X;
Xqa = Xqa./(ones(4,1)*X(4,:));
Pqa = zeros(size(P));
C = zeros(4, nFrame);
for i = 1:nFrame
	Pqa(:,:,i) = P(:,:,i)/Hqa;
	Pqa(:,:,i) = Pqa(:,:,i) / det(Pqa(:,1:3,i))^(1/3);
	C(:,i) = null(Pqa(:,:,i));
	C(:,i) = C(:,i)/C(4,i);
end
l = ones(4,1); u = l;
for i = 1:3
	cvx_begin quiet
		variable v(4);
		minimize(v(i));
		subject to
			v(4) == 1;
			Xqa'*v >= 0;
			C'*v >= 0;
	cvx_end
	cvx_clear;
	l(i) = v(i);
	cvx_begin quiet
		variable v(4);
		maximize(v(i));
		subject to
			v(4) == 1;
			Xqa'*v >= 0;
			C'*v >= 0;
	cvx_end
	cvx_clear;
	u(i) = v(i);
end
end