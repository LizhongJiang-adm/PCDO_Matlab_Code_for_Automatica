function dx = OdeSnstvSqntUB_CstCtrl_PFBF(t,x,Dec,Index,N,CurGrd)
    % 参数定义
kL=0.006;mu=0.11;
Yxs=0.47;
theta=0.004;Yp=1.2;kI=0.1;Mx=0.029;Sl=400;kXP=0.01;
kP=0.0001;

    dx = zeros(4 + N*1*4,1);

	Dec = reshape(Dec,1,1);

	u(1,1) = Dec(1,1); 	du_dp(1,1) = 1;
	dudpVct(:,1) = zeros(1*N,1); dudpVct((Index-1)*1+1,1) = 1;

    % 微分方程定义
	dx(1) =(mu*x(1)*x(2))/(x(2) + kL*x(1)) - (u(1)*x(1))/x(4);
	dx(2) =(u(1)*(Sl - x(2)))/x(4) - Mx*x(1) - (theta*x(1)*x(2))/(Yp*(kP + x(2) + x(2)^2/kI)) - (mu*x(1)*x(2))/(Yxs*(x(2) + kL*x(1)));
	dx(3) =(theta*x(1)*x(2))/(kP + x(2) + x(2)^2/kI) - (u(1)*x(3))/x(4) - kXP*x(3);
	dx(4) =u(1);
	df_dx = zeros(4,4);
	df_dx(1,1) = (mu*x(2))/(x(2) + kL*x(1)) - u(1)/x(4) - (kL*mu*x(1)*x(2))/(x(2) + kL*x(1))^2;
	df_dx(1,2) = (mu*x(1))/(x(2) + kL*x(1)) - (mu*x(1)*x(2))/(x(2) + kL*x(1))^2;
	df_dx(1,3) = 0;
	df_dx(1,4) = (u(1)*x(1))/x(4)^2;
	df_dx(2,1) = (kL*mu*x(1)*x(2))/(Yxs*(x(2) + kL*x(1))^2) - (mu*x(2))/(Yxs*(x(2) + kL*x(1))) - (theta*x(2))/(Yp*(kP + x(2) + x(2)^2/kI)) - Mx;
	df_dx(2,2) = (mu*x(1)*x(2))/(Yxs*(x(2) + kL*x(1))^2) - (mu*x(1))/(Yxs*(x(2) + kL*x(1))) - (theta*x(1))/(Yp*(kP + x(2) + x(2)^2/kI)) - u(1)/x(4) + (theta*x(1)*x(2)*((2*x(2))/kI + 1))/(Yp*(kP + x(2) + x(2)^2/kI)^2);
	df_dx(2,3) = 0;
	df_dx(2,4) = -(u(1)*(Sl - x(2)))/x(4)^2;
	df_dx(3,1) = (theta*x(2))/(kP + x(2) + x(2)^2/kI);
	df_dx(3,2) = (theta*x(1))/(kP + x(2) + x(2)^2/kI) - (theta*x(1)*x(2)*((2*x(2))/kI + 1))/(kP + x(2) + x(2)^2/kI)^2;
	df_dx(3,3) = - kXP - u(1)/x(4);
	df_dx(3,4) = (u(1)*x(3))/x(4)^2;
	df_dx(4,1) = 0;
	df_dx(4,2) = 0;
	df_dx(4,3) = 0;
	df_dx(4,4) = 0;
	df_dC = zeros(4,1);
	df_dC(1,1) = -x(1)/x(4);
	df_dC(2,1) = (Sl - x(2))/x(4);
	df_dC(3,1) = -x(3)/x(4);
	df_dC(4,1) = 1;

    for i = 1:N
		Idx0 = 4 + (i-1)*4*1;
        if i == Index
			dx(Idx0+1) = df_dx(1,1)*x(Idx0+1) +  df_dx(1,2)*x(Idx0+2) +  df_dx(1,4)*x(Idx0+4) + df_dC(1,1) *du_dp(1,1) ;
			dx(Idx0+2) = df_dx(2,1)*x(Idx0+1) +  df_dx(2,2)*x(Idx0+2) +  df_dx(2,4)*x(Idx0+4) + df_dC(2,1) *du_dp(1,1) ;
			dx(Idx0+3) = df_dx(3,1)*x(Idx0+1) +  df_dx(3,2)*x(Idx0+2) +  df_dx(3,3)*x(Idx0+3) +  df_dx(3,4)*x(Idx0+4) + df_dC(3,1) *du_dp(1,1) ;
			dx(Idx0+4) = 0 + df_dC(4,1) *du_dp(1,1) ;

        else
			dx(Idx0+1) = df_dx(1,1)*x(Idx0+1) +  df_dx(1,2)*x(Idx0+2) +  df_dx(1,4)*x(Idx0+4) ;
			dx(Idx0+2) = df_dx(2,1)*x(Idx0+1) +  df_dx(2,2)*x(Idx0+2) +  df_dx(2,4)*x(Idx0+4) ;
			dx(Idx0+3) = df_dx(3,1)*x(Idx0+1) +  df_dx(3,2)*x(Idx0+2) +  df_dx(3,3)*x(Idx0+3) +  df_dx(3,4)*x(Idx0+4) ;
			dx(Idx0+4) = 0 ;

        end
    end

	Idx0 = 4 + N*4*1 + 0*(1+1*N);
    		dx(Idx0 + 1 ) = Smoothing_Max_with_0(dx(2),0.001);
    		dDotgdp = dx(4+2:4:4+4*N);
    		dx(Idx0 + (2:N+1) ) = Smoothing_Max_with_0_der(dx(2),0.001)*dDotgdp;

end