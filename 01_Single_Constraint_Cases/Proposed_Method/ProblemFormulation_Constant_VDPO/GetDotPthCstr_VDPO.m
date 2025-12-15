function dotg = GetDotPthCstr_VDPO(t,x,u)
% 自动生成的路径约束导数值



% 计算路径约束的导数
dotg = x(:,2) - u(:,1) + x(:,1).*(x(:,2).^2 - 1);

end