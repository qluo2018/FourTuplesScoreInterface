% find solvent accessibility changes by residue number and chain ID
function saChange = findPSAbyResNum(number, chain, Joypsa)
saChange = 0;
for i = 1 : size(Joypsa,2)
    if Joypsa(i).resnum == number ...
            && Joypsa(i).reschain == chain 
        saChange = Joypsa(i).sachange;
        break;
    end
end