function index = MapResiduesToIndex20(residue, ID)
index = -1;
for i = 1 : size(ID,2)
    if strcmp(residue,ID{i})
        index = i;
        break;
    end
end
