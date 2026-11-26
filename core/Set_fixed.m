%% Overwrites the fields of p with the fixed hyperparameters
function p = Set_fixed(p, fixed)
    names = fieldnames(fixed);
    for i = 1:numel(names)
        p.(names{i}) = fixed.(names{i});
    end
end
