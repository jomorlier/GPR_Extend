function new_kern = pseudo_kern_trick(OM)
%Performs monomial higher-order projection to the second degree, results in
%products of all cols with each other and squares
    new_kern = [];
    for i = 1:length(OM(1,:))
        for j = i:length(OM(1,:))
            new_kern =  [new_kern OM(:,i).*OM(:,j)];
        end
    end
end