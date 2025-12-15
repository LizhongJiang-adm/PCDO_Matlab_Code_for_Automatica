function [Cieq, Ceq, G_Cieq, G_Ceq] = Constraints_Dis_JCB(u,CvpGrd,x0)
% 自动生成的多重打靶法约束文件，可以包含预设的离散等式约束和离散不等式约束

N = numel(CvpGrd)-1;
NStat = numel(x0);
NPmtr = 1; NCtrl = 1;
NCPm = NPmtr*NCtrl;
NDec1 = NCPm;  % 一个控制区间上的决策变量数量
NDec = N*NDec1;       % 总决策变量数
u = reshape(u,numel(u),1);
x0 = reshape(x0,1,numel(x0));
TranGrd = CvpGrd; % Cvp分段时间点
NTrGd = numel(TranGrd)-1;

DesMat = reshape(u,NDec1,N); % 第一个区间初值也是决策变量

XendAug = [x0, zeros(1,NStat*NDec)];
GraXend = zeros(NDec,NStat);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

% 定义离散路径约束
CieqDisPath1 = zeros(1,N+1);
G_CieqDisPath1 = zeros(NDec,N+1);

% 对初始时刻进行约束
StatTD = XendAug; % 插值求得系统状态在TD的值
GraStatTD = reshape(XendAug(NStat+1:NStat+NStat*NDec),NStat,NDec).';
CieqDisPath1(1,1) = StatTD(2) - 8*(CvpGrd(1) - 1/2)^2 + 1/2;
G_CieqDisPath1(:,1) = GraStatTD(:,2);

for i = 1:NTrGd

	C_t0 = TranGrd(i);  C_tf = TranGrd(i+1);
	for ii=1:N % 确定当前控制阶段，取出对应决策变量 
		if(CvpGrd(ii)-C_t0<1e-10) && (C_tf-CvpGrd(ii+1)<1e-10)
		CurIdx = ii;   CurDec = DesMat(:,CurIdx); 
		end
	end

	StatAug0 = XendAug; % 拼接系统状态

	% 对系统进行积分
	[t, x] = ode45(@(t,x)OdeSnstvSqnt_CstCtrl_JCB(t,x,CurDec,CurIdx,N,[C_t0,C_tf]),...
	[C_t0,C_tf], StatAug0, opts);

	XendAug = x(end,:);

	GraXend = reshape(XendAug(NStat+1:NStat+NStat*NDec),NStat,NDec).';

	% 取出离散约束

	StatTD = XendAug; % 插值求得系统状态在TD的值
	GraStatTD = GraXend;
	CieqDisPath1(1,i+1) = StatTD(2) - 8*(C_tf - 1/2)^2 + 1/2;
	G_CieqDisPath1(:,i+1) = GraStatTD(:,2);

	% 取出离散路径约束
	% for j=1:numel(T_DisPath1)
	% 	TD = T_DisPath1(j);
	% 	if ( (TD>=TranGrd(i)) && ( (i==N) || (TD<TranGrd(i+1))))==1 % 找出区间中的TD
	% 		Idx0a = 0;
	% 		StatTD = interp1(t,x,TD,'PCHIP'); % 插值求得系统状态在TD的值
	% 		GraStatTD = reshape(StatTD(NStat+1:NStat+NStat*NDec),NStat,NDec).';
	% 		CieqDisPath1(1,j) = StatTD(2) - 8*(TD - 1/2)^2 + 1/2;
	% 		G_CieqDisPath1(:,j) = GraStatTD(:,2);
	% 	end
	% end

end

Ceq = [];
G_Ceq = [];
Cieq = [ ];
G_Cieq = [ ];
Cieq = [Cieq, CieqDisPath1];
G_Cieq = [G_Cieq, G_CieqDisPath1];

end