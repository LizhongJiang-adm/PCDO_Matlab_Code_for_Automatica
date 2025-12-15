function dx = OdeSnstvSqnt_CstCtrl_RAYL(t,x,Dec,Index,N,CurGrd)
    % 参数定义


    dx = zeros(3 + N*1*3,1);

	Dec = reshape(Dec,1,1);

	u(1,1) = Dec(1,1); 	du_dp(1,1) = 1;
	dudpVct(:,1) = zeros(1*N,1); dudpVct((Index-1)*1+1,1) = 1;

    % 微分方程定义
	dx(1) =x(2);
	dx(2) =4*u(1) - x(1) - x(2)*((7*x(2)^2)/50 - 7/5);
	dx(3) =u(1)^2 + x(1)^2;



    for i = 1:N
		Idx0 = 3 + (i-1)*3*1;
        if i == Index
			dx(Idx0+1) = (1)*x(Idx0+2) + 0 ;
			dx(Idx0+2) = (-1)*x(Idx0+1) +  (7/5 - (21*x(2)^2)/50)*x(Idx0+2) + (4) *du_dp(1,1) ;
			dx(Idx0+3) = (2*x(1))*x(Idx0+1) + (2*u(1)) *du_dp(1,1) ;

        else
			dx(Idx0+1) = (1)*x(Idx0+2) ;
			dx(Idx0+2) = (-1)*x(Idx0+1) +  (7/5 - (21*x(2)^2)/50)*x(Idx0+2) ;
			dx(Idx0+3) = (2*x(1))*x(Idx0+1) ;

        end
end