function index = MapResiduesToIndex(residue, ID)
AAtype = [6	5 3	2 5	3 1	6 2	2 7	2 4	1 1	5 4	4 2	1];
index = -1;
for i = 1 : size(ID,2)
    if strcmp(residue,ID{i})
        index = AAtype(i);
        break;
    end
end
