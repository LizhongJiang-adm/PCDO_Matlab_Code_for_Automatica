



max_theta = 0.25; max_speed = 1.5; max_phi = 0.714;

% 计算状态轨迹与最大值
[X_sol, U_sol, ~, T_sol] = strategy.unpack_solution(w_opt);
t_dense = linspace(0, T_sol, 20000);
[x_dense, u_dense] = strategy.get_dense_trajectory(w_opt, t_dense);
[MaxValTraj, IdxMax]= max(x_dense(3:5,:),[],2);
MaxTimeTraj = t_dense(IdxMax);

% 绘图参数
x_Range = [max_theta,max_speed,max_phi];
LgndTrajectory = {'$x_3(t)$' '$x_4(t)$' '$x_5(t)$'};
LgndRange = {'State Boundaries' 'State Boundaries' 'State Boundaries'};
LgndMax = {'$\max_{t\in\Gamma}x_3(t)$' '$\max_{t\in\Gamma}x_4(t)$' '$\max_{t\in\Gamma}x_5(t)$'};
LgndPos = [0.124121325327156,0.036025635303595,0.769688198482368,0.111153841972351];

TMv = [1e-4, 5e-5, 5e-5];     GMv = [5e-6, 1e-4, 5e-5];

SubPicPst = cell(1,3);
SubPicPst{1} = [0.436404764007767,0.319047619047622,0.295738093135089,0.169270320681385];
SubPicPst{2} = [0.436404764007767,0.319047619047622,0.295738093135089,0.169270320681385];
SubPicPst{3} = [0.450690478293483,0.523809523809524,0.295738093135089,0.169270320681384];

ArrPst = {[0.2946,0.4054;0.669,0.5333]    [0.3071,0.4018;0.6833,0.5381]    [0.2464,0.3786;0.6762,0.6238]};

FileName = {"StateTrajectory_Parking_x3.pdf", "StateTrajectory_Parking_x4.pdf", "StateTrajectory_Parking_x5.pdf"};

for i = 1:3
    figure
    hold on
    plot(t_dense,x_dense(i+2,:),'LineWidth',2,DisplayName=LgndTrajectory{i});
    yline(x_Range(i),'--',[],"Color",'m','LineWidth',2,DisplayName=LgndRange{i});
    yline(-x_Range(i),'--',[],"Color",'m','LineWidth',2,HandleVisibility='off');
    plot(MaxTimeTraj(i), MaxValTraj(i), 'p',...
        'MarkerFaceColor','b','MarkerEdgeColor','b','markersize', 8, DisplayName=LgndMax{i});
    
    annotation('textarrow',ArrPst{i}(1,:),ArrPst{i}(2,:));

    YLimMid = 0;
    Gap = 2*x_Range(i);
    ylim([YLimMid-0.7*Gap, YLimMid+0.9*Gap])
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 15);
    legend('Interpreter','latex','FontSize',13);

    box on
    hold off


    axes('Position',SubPicPst{i});
    box on
    set(groot,'defaultAxesLineStyleOrder','remove','defaultAxesColorOrder','remove')
    
    hold on
    plot(t_dense,x_dense(i+2,:),'LineWidth',2,DisplayName=LgndTrajectory{i});
    yline(x_Range(i),'--',[],"Color",'m','LineWidth',2,DisplayName=LgndRange{i});
    yline(-x_Range(i),'--',[],"Color",'m','LineWidth',2,HandleVisibility='off');
    plot(MaxTimeTraj(i), MaxValTraj(i), 'p','MarkerFaceColor','b','MarkerEdgeColor','b','markersize', 8, DisplayName=LgndMax{i});    
    hold off

    YLimMid = 0.5*(MaxValTraj(i) + x_Range(i));
    Gap = x_Range(i) - MaxValTraj(i);
    xlim([MaxTimeTraj(i)-TMv(i), MaxTimeTraj(i)+TMv(i)]);
    % ylim([YLimMid-GMv(i), YLimMid+GMv(i)])
    ylim([YLimMid-0.8*Gap, YLimMid+0.8*Gap])

    if saveResult
        exportgraphics(gcf, FileName{i}, 'Resolution', 600, 'ContentType', 'vector');
    end
end


% 
% 
% ArrPst = {[0.28,0.5586;0.6262,0.5215]    [0.3071,0.4018;0.6833,0.5381]    [0.1843,0.2257;0.7108,0.5246]};
% ColorList = ["#0072BD","#D95319","#EDB120"];
% 
% 
% figure(Position=[781,431.4,560,520])
% axes('Position',[0.122857142857143,0.218461538461538,0.775,0.748076923076924]);
% box on
% hold on
% for i = 1:3
%     plot(t_dense,x_dense(i+2,:),"Color",ColorList(i),'LineWidth',1.2,DisplayName=LgndTrajectory{i});
%     yline(x_Range(i),'--',[],"Color",ColorList(i),'LineWidth',1.5,DisplayName=LgndRange{i});
%     yline(-x_Range(i),'--',[],"Color",ColorList(i),'LineWidth',1.5,HandleVisibility='off');
%     plot(MaxTimeTraj(i), MaxValTraj(i), 'p',...
%         'MarkerFaceColor',ColorList(i),'MarkerEdgeColor',ColorList(i),'markersize', 8, DisplayName=LgndMax{i});
%     % annotation('textarrow',ArrPst{i}(1,:),ArrPst{i}(2,:));
% end
% legend('Interpreter','latex','FontSize',11,NumColumns=3,Position=LgndPos);
% ylim([-1.8,1.8])
% 
% 
% 
% 
% 
% for i = 1:3
%     axes('Position',SubPicPst{i});
%     box on
%     set(groot,'defaultAxesLineStyleOrder','remove','defaultAxesColorOrder','remove')
% 
%     hold on
%     plot(t_dense,x_dense(i+2,:),"Color",ColorList(i),'LineWidth',1.2,DisplayName=LgndTrajectory{i});
%     yline(x_Range(i),'--',[],"Color",ColorList(i),'LineWidth',1.5,DisplayName=LgndRange{i});
%     yline(-x_Range(i),'--',[],"Color",ColorList(i),'LineWidth',1.5,HandleVisibility='off');
%     plot(MaxTimeTraj(i), MaxValTraj(i), 'p','MarkerFaceColor',ColorList(i),'MarkerEdgeColor',ColorList(i),'markersize', 8, DisplayName=LgndMax{i});
% 
% 
%     hold off
% 
%     YLimMid = 0.5*(MaxValTraj(i) + x_Range(i));
%     Gap = x_Range(i) - MaxValTraj(i);
%     xlim([MaxTimeTraj(i)-TMv(i), MaxTimeTraj(i)+TMv(i)]);
%     % ylim([YLimMid-GMv(i), YLimMid+GMv(i)])
%     ylim([YLimMid-0.8*Gap, YLimMid+0.8*Gap])
% end
