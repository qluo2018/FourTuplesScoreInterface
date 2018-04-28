function joy = findjoypsaforresidue(joypsa, id, residue)
joy = [];
for i = 1 : size(joypsa,2)
    if residue.number(id) == joypsa(i).resnum && strcmpi(residue.chain(id), joypsa(i).reschain)
        joy = joypsa(i);
    end
end