function dotg = GetDotPthCstr_PFBF(t,x,u)
% 自动生成的路径约束导数值

kL=0.006;mu=0.11;
Yxs=0.47;
theta=0.004;Yp=1.2;kI=0.1;Mx=0.029;Sl=400;kXP=0.01;
kP=0.0001;

% 计算路径约束的导数
dotg = (u(:,1).*(Sl - x(:,2)))./x(:,4) - Mx.*x(:,1) - (theta.*x(:,1).*x(:,2))./(Yp.*(kP + x(:,2) + x(:,2).^2./kI)) - (mu.*x(:,1).*x(:,2))./(Yxs.*(x(:,2) + kL.*x(:,1)));

end