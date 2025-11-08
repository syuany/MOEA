function PlotVals(pop)
    figure(1);
    Values = [pop.Value];
    plot(Values(1, :), Values(2, :), 'x');
    xlabel('1^{st} Objective');
    ylabel('2^{st} Objective');
    grid on;
end