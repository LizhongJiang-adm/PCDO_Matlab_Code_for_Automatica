function dotg = GetDotPthCstr_JCB(t,x,u)
% 自动生成的路径约束导数值



% 计算路径约束的导数
dotg = u(:,1) - 16.*t - x(:,2) + 8;

end