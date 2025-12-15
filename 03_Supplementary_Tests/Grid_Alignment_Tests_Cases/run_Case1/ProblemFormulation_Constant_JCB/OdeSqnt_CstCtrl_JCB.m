function dx = OdeSqnt_CstCtrl_JCB(t,x,Dec,Index,N,CurGrd)
    % 参数定义


    dx = zeros(2,1);

	Dec = reshape(Dec,1,1);

	u(1,1) = Dec(1,1); 	du_dp(1,1) = 1;
	dudpVct(:,1) = zeros(1*N,1); dudpVct((Index-1)*1+1,1) = 1;

    % 微分方程定义
	dx(1) =x(2);
	dx(2) =u(1) - x(2);



end