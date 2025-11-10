function PlotCosts(pop)
    figure(1);
    Costs = [pop.Cost];
    plot(Costs(1, :), Costs(2, :), 'x');
    xlabel('1^{st} Objective');
    ylabel('2^{st} Objective');
    grid on;
end