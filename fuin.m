function x = fuin(x, l, u)
cvx_begin
	subject to
		l <= x <= u;
cvx_end
end