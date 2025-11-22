function Population = GDE3_EnvironmentalSelection(Population, Offspring, N)
    PopCost=[Population.Cost];
    OffCost=[Offspring.Cost];
    updated=all(PopCost>OffCost, 2);
    selected=any(PopCost>OffCost, 2);
    Population(updated)=Offspring(updated);
    Population=[Population; Offspring(selected)];
    Population=EnvironmentalSelection(Population, N);
end