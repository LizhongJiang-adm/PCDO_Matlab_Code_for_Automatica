function [objval,Gobj] = ObjtvSqnt_CstCtrl_PFBF(u,CvpGrd,x0)
% 自动生成的多重打靶法目标函数文件，不含积分目标版本

N = numel(CvpGrd)-1;
NStat = numel(x0);
NPmtr = 1; NCtrl = 1;
NCPm = NPmtr*NCtrl;
NDec1 = NCPm;  % 一个控制区间上的决策变量数量
NDec = numel(u);       % 总决策变量数

DesMat = reshape(u,NDec1,N); % 第一个区间初值也是决策变量
StatAug0 = [x0,zeros(1,NStat*NDec)];

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

for i = 1:N
	CurIdx = i; CurDec = DesMat(:,i);
	C_t0 = CvpGrd(CurIdx);  C_tf = CvpGrd(CurIdx+1);
	[t, x] = ode45(@(t,x)OdeSnstvSqntObj_CstCtrl_PFBF(t,x,CurDec,CurIdx,N,[C_t0,C_tf]),...
	[C_t0,C_tf], StatAug0, opts);

	StatAug0 = x(end,:);
end

objval = 0;
Gobj = zeros(NDec,1);

XendAug = x(end,:);
GraXend = reshape(XendAug(NStat+1:NStat+NStat*NDec),NStat,NDec).';

objval = objval + -XendAug(3);
Gobj = Gobj + -GraXend(:,3);

end