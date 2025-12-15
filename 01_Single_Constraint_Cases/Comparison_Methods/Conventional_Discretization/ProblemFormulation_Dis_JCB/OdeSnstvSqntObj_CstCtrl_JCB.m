function dx = OdeSnstvSqntObj_CstCtrl_JCB(t,x,Dec,Index,N,CurGrd)
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

	Idx0 = 2 + N*2;
	dx(Idx0+1) = u(1)^2/200 + x(1)^2 + x(2)^2;
	dx(Idx0+2:Idx0+1*N+1) = (dudpVct(:,1)*u(1))/100 + 2*dxdpVct(:,1)*x(1) + 2*dxdpVct(:,2)*x(2);
end