function [l_mat,u_mat] = prop_interval(l, u)
l_k = [l(1:3)'; 0 l(4:5)'; 0 0 1];
u_k = [u(1:3)'; 0 u(4:5)'; 0 0 1];
assert(sign(l_k(1,1)) == sign(u_k(1,1)));
assert(sign(l_k(2,2)) == sign(u_k(2,2)));
assert(sign(l_k(1,2)) ~= sign(u_k(1,2)));
assert(sign(l_k(1,3)) ~= sign(u_k(1,3)));
assert(sign(l_k(2,3)) ~= sign(u_k(2,3)));
l_11 = l_k(1,1) * l_k(1,1);
l_12 = l_k(1,2) * u_k(2,2) + min(l_k(1,3) * u_k(2,3), u_k(1,3) * l_k(2,3));
l_13 = l_k(1,3);
l_21 = l_12;
l_22 = l_k(2,2) * l_k(2,2);
l_23 = l_k(2,3);
l_31 = l_13;
l_32 = l_23;
l_33 = 1;

l_mat = [[l_11, l_12, l_13];[l_21,l_22,l_23];[l_31,l_32,l_33]];

u_11 = u_k(1,1) * u_k(1,1) + max(l_k(1,2) * l_k(1,2), u_k(1,2) * u_k(1,2)) + max(l_k(1,3) * l_k(1,3), u_k(1,3)*u_k(1,3));
u_12 = u_k(1,2) * u_k(2,2) + max(l_k(1,3) * l_k(2,3), u_k(1,3) * u_k(2,3));
u_13 = u_k(1,3);
u_21 = u_12;
u_22 = u_k(2,2) * u_k(2,2) + max(l_k(2,3) * l_k(2,3), u_k(2,3) * u_k(2,3));
u_23 = u_k(2,3);
u_31 = u_13;
u_32 = u_23;
u_33 = 1;

u_mat = [[u_11,u_12,u_13];[u_21,u_22,u_23];[u_31,u_32,u_33]];
end