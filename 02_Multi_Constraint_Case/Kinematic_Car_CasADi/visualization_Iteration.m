

figure 
hold on

NumIteration = numel(stats.history.iterationPoint);

for i = 1:NumIteration
    [PathCnstrMaxVal, PathCnstrMaxTime_nmlz, ~, ~] = strategy.getMaxPConIntvls(stats.history.iterationPoint{i});
    for j = 1:numel(problem.path_constraints)
        MaxPathIteration(j,i) = max(PathCnstrMaxVal{j});
    end

    FminIteration(i) = stats.history.objVal{i};
end



% 左侧纵坐标轴
yyaxis left;
% plot(1:NumIteration, MaxPathIteration([1 3 5],:), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);

MaxPathIteration = max(MaxPathIteration);
plot(1:NumIteration, MaxPathIteration, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);

ylabel('Maximum Value Across All Path Constraints');  % 左纵轴标签
grid on;
plot(1:NumIteration,zeros(1,NumIteration),'--','LineWidth',1.5);


% RangeMaxPath = [min(MaxPathIteration),0];
range_y = ylim;
min_y = min(MaxPathIteration);
max_y = 0;
height_y = max_y - min_y;
padding_down = height_y * 0.7;
padding_up = height_y * 0.6;
% 设置新的 Y 轴范围
ylim([min_y - padding_down, max_y + padding_up]);
% ylim(PathLim)


% 右侧纵坐标轴
yyaxis right;
plot(1:NumIteration, FminIteration, 's-', 'LineWidth', 1.5, 'MarkerSize', 6);
ylabel('Objective Value');      % 右纵轴标签

% 公共横坐标轴标签和标题
xlabel('Iteration Count');
xticks(1:NumIteration);

legend({"$\max_{t\in\Gamma,j\in \{1,...,6\}}g_j(u^k,t)$","$\max_{t\in\Gamma}g(u^k,t)=0$","$h(u^k)$"},'Interpreter','latex','FontSize',11);


hold off
box on

if saveResult
    exportgraphics(gcf, "IterationTrend.pdf", 'Resolution', 600, 'ContentType', 'vector');
end