function dx = OdeSnstvSqntUB_CstCtrl_JCB(t,x,Dec,Index,N,CurGrd)
    % 参数定义


    dx = zeros(2 + N*1*2,1);

	Dec = reshape(Dec,1,1);

	u(1,1) = Dec(1,1); 	du_dp(1,1) = 1;
	dudpVct(:,1) = zeros(1*N,1); dudpVct((Index-1)*1+1,1) = 1;

    % 微分方程定义
	dx(1) =x(2);
	dx(2) =u(1) - x(2);



    for i = 1:N
		Idx0 = 2 + (i-1)*2*1;
        if i == Index
			dx(Idx0+1) = (1)*x(Idx0+2) + 0 ;
			dx(Idx0+2) = (-1)*x(Idx0+2) + (1) *du_dp(1,1) ;

        else
			dx(Idx0+1) = (1)*x(Idx0+2) ;
			dx(Idx0+2) = (-1)*x(Idx0+2) ;

        end
    end

	dxdpVct = reshape(x(2+1:2+N*2*1),2,1*N).';
	Idx0 = 2 + N*2*1 + 0*(1+1*N);
    		dx(Idx0 + 1 ) = Smoothing_Max_with_0(u(1) - 16*t - x(2) + 8,0.001);
    		dDotgdp = (-1)*dxdpVct(:,2) +  (1)*dudpVct(:,1) ;
    		dx(Idx0 + (2:1*N+1) ) = Smoothing_Max_with_0_der(u(1) - 16*t - x(2) + 8,0.001)*dDotgdp;

end