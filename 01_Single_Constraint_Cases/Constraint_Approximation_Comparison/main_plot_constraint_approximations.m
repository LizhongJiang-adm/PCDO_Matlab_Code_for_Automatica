

clear

addpath("PicErrFun");
%% --- 1. Define the Optimal Control Problem ---
fprintf('Step 1: Defining the PFBF optimal control problem...\n');
problem = define_PFBF_for_visualization();

load("TestPoint.mat");

global LineColr LineWdt
LineWdt = 1.5;
LineColr{1} = 'k';%[0.93,0.69,0.13];
LineColr{2} = 'r';%[0.64,0.08,0.18];
LineColr{3} = [0.47,0.67,0.19];
LineColr{4} = [0.49,0.18,0.56];


N = problem.control.N;                         %% 控制向量的网格数
T = problem.time.TF;                         %% 时间范围
CstrIntvl = problem.pathConstraints.constraintGrid;
x0 = problem.state.x0;
CvpGrd = problem.control.grid;


Pic = figure;

%% 设定绘图参数


h1=axes('Position',[0.13,0.24,0.775,0.685]);
ylim([-0.5,0.5]);
xlabel('$t$','Interpreter','latex');
ylabel('$g(u,t)$','Interpreter','latex');
box on
hold on

% 绘制路径约束曲线
[Gmax, Tmax,st_x_g,st_t,ctrl_t] = PathCnstrProfileV3(DecVar,problem); 
g_val = st_x_g(:,end);
plot(st_t,g_val,'b-.','LineWidth',2,DisplayName=sprintf('$g(u,t)$'));
plot(st_t,st_x_g(:,1)-st_x_g(:,1),'m--', DisplayName='$g(u,t)=0$');


% 绘制约束区间网格——灰色虚线
for i = 1:numel(CstrIntvl)
    Xpnt = CstrIntvl(i)*[1,1];
    CstrGrd = plot(Xpnt,ylim,'--','Color',[0.7, 0.7, 0.7],'HandleVisibility','off', 'LineWidth',0.8);
end

[Cieq_UB, ~, ~, ~] = CnstrUBSqnt_CstCtrl_PFBF(DecVar,CvpGrd,x0,CstrIntvl);
plot_horizontal_lines(CstrIntvl,Cieq_UB,'-',Color=LineColr{1},LineWidth=LineWdt,DisplayName='$S(u,T_i)$');

ResPmtr = zeros(1,N);
[Cieq_Poly, ~, ~, ~] = Constraints_Poly_PFBF(DecVar,problem, CstrIntvl, ResPmtr);
plot_horizontal_lines(CstrIntvl,Cieq_Poly(N+2:end),'-',Color=LineColr{2},LineWidth=LineWdt,DisplayName='$g^{\max}_{\mathrm{Poly}}(u,t)$');

[Cieq_Intvl, ~, ~, ~] = Constraints_Intvl_PFBF(DecVar,CstrIntvl,problem);
plot_horizontal_lines(CstrIntvl,Cieq_Intvl,'-',Color=LineColr{3},LineWidth=LineWdt,DisplayName='$g^{\max}_\mathrm{Intvl}(u,T_i)$');

t_Rlt = [];
myalpha = problem.algorithm_params.myalpha;
for i = 1:N
    Idx = find(CvpGrd(i)<=st_t & st_t<CvpGrd(i+1));
    t_Rlt = st_t(Idx)-CvpGrd(i); % [t_Rlt;st_t(Idx)-CvpGrd(i)];
    [MaBB,~] = max(g_val(Idx,1) + myalpha*t_Rlt.^2);
    plot([CvpGrd(i),CvpGrd(i+1)], [MaBB,MaBB],'-',Color=LineColr{4},LineWidth=LineWdt,HandleVisibility='off');
end
plot([CvpGrd(i),CvpGrd(i+1)], [MaBB,MaBB],'-',Color=LineColr{4},LineWidth=LineWdt,DisplayName='$g^{\max}_{\mathrm{\alpha BB}}(u,T_i)$');

LgndGmax = ['$\max\limits_{t\in T_i} g(u,t)$'];
plot(Tmax, Gmax, 'bp','MarkerFaceColor','b','markersize', 8, DisplayName=LgndGmax);

HLgd = legend('Interpreter','latex','FontSize',12,NumColumns=4);
set(HLgd,'Position',[0.146347216652202,0.026403830389211,0.732543318547081,0.105223154475737])




SubPicPst = [0.241428571428569,0.308253962119305,0.269274193548385,0.189841275975937];
h2=axes('Position',SubPicPst);
hold on

set(groot,'defaultAxesLineStyleOrder','remove','defaultAxesColorOrder','remove')
plot(st_t,g_val,'b-.','LineWidth',2,DisplayName=sprintf('$g(u,t)$'));
plot(st_t,st_x_g(:,1)-st_x_g(:,1),'m--', DisplayName='$g(u,t)=0$');
plot_horizontal_lines(CstrIntvl,Cieq_UB,'-',Color=LineColr{1},LineWidth=LineWdt,DisplayName='$S(u,T_i)$');
plot_horizontal_lines(CstrIntvl,Cieq_Poly(N+2:end),'-',Color=LineColr{2},LineWidth=LineWdt,DisplayName='$g^{\max}_{\mathrm{Poly}}(u,t)$');
plot(Tmax, Gmax, 'bp','MarkerFaceColor','b','markersize', 8, DisplayName=LgndGmax);

ArrPst = [0.5339,0.4393;0.5571,0.5119];
annotation('textarrow',ArrPst(1,:),ArrPst(2,:));

PicLim = [Tmax-5e-4, Tmax+5e-4; -0.1e-3,0.3e-3];
set(h2,'xlim',PicLim(1,:),'ylim',PicLim(2,:));
box on

exportgraphics(gcf, "constraint_approximation_comparison_PFBF.pdf", 'Resolution', 600, 'ContentType', 'vector');
