function [Cieq, Ceq, G_Cieq, G_Ceq] = CnstrUBSqnt_CstCtrl_PFBF(u,CvpGrd,x0,CstrIntvl)
% 自动生成的多重打靶法约束文件，可以包含预设的离散等式约束和离散不等式约束

N = numel(CvpGrd)-1;
NStat = numel(x0);
NPmtr = 1; NCtrl = 1;
NCPm = NPmtr*NCtrl;
NDec1 = NCPm;  % 一个控制区间上的决策变量数量
NDec = N*NDec1;       % 总决策变量数
u = reshape(u,numel(u),1);
x0 = reshape(x0,1,numel(x0));
TranGrd = unique([CvpGrd,CstrIntvl]); % 在Cvp节点上切换控制信号，在CstrIntvl节点上取出约束，并更新积分约束初值
NTrGd = numel(TranGrd)-1;

DesMat = reshape(u,NDec1,N); % 第一个区间初值也是决策变量

XendAug = [x0, zeros(1,NStat*NDec)];
GraXend = zeros(NDec,NStat);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

% 定义离散约束

% 定义上界约束
CstrIntvl1 = CstrIntvl;
Cieq_UB1 = zeros(1,numel(CstrIntvl1)-1);
G_Cieq_UB1 = zeros(NDec,numel(CstrIntvl1)-1);

for i = 1:NTrGd

	C_t0 = TranGrd(i);  C_tf = TranGrd(i+1);
	for ii=1:N % 确定当前控制阶段，取出对应决策变量 
		if(CvpGrd(ii)-C_t0<1e-10) && (C_tf-CvpGrd(ii+1)<1e-10)
			CurIdx = ii;   CurDec = DesMat(:,CurIdx); 
			CurCVPGrd = [CvpGrd(ii),CvpGrd(ii+1)];
		end
	end

	StatAug0 = XendAug; % 拼接系统状态

	IdxC_t0 = find(C_t0==CstrIntvl1,1);
	if ~isempty(IdxC_t0) % 如果当前阶段起始点是某个约束区间的网格，更新积分约束初值
		Idx0a = NStat+NStat*NDec + 0*(1+NDec1);
		StatAug0(Idx0a+1) = StatAug0(2) - 1/2;

		GraStatAug0 = reshape(StatAug0(NStat+1:NStat+NStat*NDec),NStat,NDec).';
		StatAug0(Idx0a+2:Idx0a+1+NDec) = GraStatAug0(:,2);
	end

	% 对系统进行积分
	[t, x] = ode45(@(t,x)OdeSnstvSqntUB_CstCtrl_PFBF(t,x,CurDec,CurIdx,N,CurCVPGrd),...
	[C_t0,C_tf], StatAug0, opts);

	XendAug = x(end,:);

	GraXend = reshape(XendAug(NStat+1:NStat+NStat*NDec),NStat,NDec).';

	% 取出离散约束

	% 取出上界约束
	IdxC_tf = find(C_tf==CstrIntvl1,1);
	if ~isempty(IdxC_tf) && IdxC_tf~=1 
	% 如果当前阶段终点是某个约束区间的网格，而且不是第一个区间左端点（此时约束没有经过积分），理论上不可能出现，则取出约束值
		Idx0a = NStat+NStat*NDec + 0*(1+NDec);
		Cieq_UB1(1,IdxC_tf-1) = XendAug(Idx0a+1);
		G_Cieq_UB1(:,IdxC_tf-1) = XendAug(Idx0a+2:Idx0a+1+NDec);
	end

end

Ceq = [];
G_Ceq = [];
Cieq = [ ];
G_Cieq = [ ];
Cieq = [Cieq, Cieq_UB1];
G_Cieq = [G_Cieq, G_Cieq_UB1];

end