% find solvent accessibility changes by residue number and chain ID
function saTotalSide = findPSATotalSidebyResNum(number, chain, Joypsa)
saTotalSide = 0;
for i = 1 : size(Joypsa,2)
    if Joypsa(i).resnum == number ...
            && strcmpi(Joypsa(i).reschain, chain) 
        saTotalSide = Joypsa(i).splitedtotalside;
        break;
    end
end