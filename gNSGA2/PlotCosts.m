function PlotCosts(pop, g)
    figure(1);  
    % gcf 是 "get current figure" 的缩写，表示获取当前的图形窗口句柄
    % 'Position' 是图形窗口的属性，用于设置窗口的位置和尺寸
    % [100 100 800 600] 是位置参数数组，包含四个值：
    % 第1个值 100：窗口左上角距离屏幕左边的距离（像素）
    % 第2个值 100：窗口左上角距离屏幕顶部的距离（像素）
    % 第3个值 800：窗口的宽度（像素）
    % 第4个值 600：窗口的高度（像素）
    set(gcf, 'Position', [100 500 600 400]);
    
    f1_exact = linspace(0, 1, 200);
    f2_exact = 1 - sqrt(f1_exact);

    % 'k-'：黑色实线
    % 'LineWidth', 2：线宽为2个像素
    % 'DisplayName', 'Exact'：在图例中显示为"Exact"
    plot(f1_exact, f2_exact, 'k-', 'LineWidth', 1, 'DisplayName', 'Exact');
    
    hold on;
    
    Costs = [pop.Cost];
    
    validIdx = all(Costs < 1000, 1);
    if any(validIdx)
        % 'rx'：绘图样式：红色的"x"标记符号
        % 'MarkerSize', 6：标记符号的大小为6
        % 'LineWidth', 1.5：线条宽度为1.5
        plot(Costs(1, validIdx), Costs(2, validIdx), 'rx', ...
                'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'NSGA-II');
    end
    
    if length(g) >= 2
        % 'bo'：绘图样式：蓝色圆圈标记
        plot(g(1), g(2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'g');
    end
    
    xlim([0 1.2]);
    ylim([0 1.5]);

    xticks(0:0.2:1.2);  % X轴步长为0.2
    yticks(0:0.2:1.5);  % Y轴步长为0.2
    
    % xlabel('1^{st} Objective', 'FontSize', 8);
    % ylabel('2^{nd} Objective', 'FontSize', 8);
    title('Pareto Front Comparison', 'FontSize', 10);
    
    legend('Location', 'northwest', 'FontSize', 8);
    
    grid on;
    box on;
    hold off;
end