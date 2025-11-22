function Population=EnvironmentalSelection(Population, N)
    [Population, F] = NDSort(Population);
    Population = CalcCrowdingDistance(Population, F);
    % Population = SortPopulation(Population);

    % Sort Based on Crowding Distance
    [~, CDSO] = sort([Population.CrowdingDistance], 'descend');
    Population = Population(CDSO);

    % Sort Based on Rank
    [~, RSO] = sort([Population.Rank]);
    Population = Population(RSO);

    if numel(Population)>N
        Population = Population(1:N);
    end
end