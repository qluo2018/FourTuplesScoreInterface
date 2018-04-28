function count_con_clique2 = CliqueCounting4forcon(sitesondomain1, s1, sitesondomain2, s2, count_con_clique, fourcliqueType, IndexSeq, structuregraph1)
%#eml
% emlmex CliqueCounting4forcon -eg {zeros(1,2000), zeros(1,1), zeros(1,2000), zeros(1,1), zeros(210,11), zeros(7,7,7,7), zeros(1,4000), zeros(2000,2000) }
% assert(all(size(sitesondomain1) == [1,2000]));
% assert(isa(sitesondomain1,'numeric'));
% assert(isa(s1, 'numeric'));
% assert(all(size(sitesondomain2) == [1,2000]));
% assert(isa(sitesondomain2,'numeric'));
% assert(isa(s2,'numeric'));
% assert(all(size(count_all_clique) == [210,11]));
% assert(isa(count_all_clique,'numeric'));
% assert(all(size(fourcliqueType) == [7,7,7,7]));
% assert(isa(fourcliqueType, 'numeric'));
% assert(all(size(IndexSeq) == [1,4000]));
% assert(isa(IndexSeq,'numeric'));
% assert(all(size(structuregraph1) == [2000,2000]));
% assert(isa(structuregraph1), 'numeric')

% pick up two sites from each domain to form the 4-clique candidate
for i = 1 : s1-1
    for j = i+1 : s1
        for p = 1 : s2-1
            for q = p+1 : s2
                fourseqindex = IndexSeq([sitesondomain1(i), sitesondomain1(j), sitesondomain2(p), sitesondomain2(q)]);
                if fourseqindex(1)> 0 && fourseqindex(2) > 0 && fourseqindex(3) > 0  && fourseqindex(4) > 0
                    cliquetyep = fourcliqueType( fourseqindex(1), fourseqindex(2), fourseqindex(3), fourseqindex(4) );
                else
                    cliquetyep = 0;
                end
                if cliquetyep > 0
                    summ = structuregraph1(sitesondomain1(i), sitesondomain1(j)) ...
                        + structuregraph1(sitesondomain1(i), sitesondomain2(p)) ...
                        + structuregraph1(sitesondomain1(i), sitesondomain2(q)) ...
                        + structuregraph1(sitesondomain1(j), sitesondomain2(p)) ...
                        + structuregraph1(sitesondomain1(j), sitesondomain2(q)) ...
                        + structuregraph1(sitesondomain2(p), sitesondomain2(q));
                    degreeofthree = getdegreeofnode(structuregraph1,sitesondomain1, sitesondomain2, i,j,p,q);
                    numofintra = countintraedge(structuregraph1,sitesondomain1, sitesondomain2, i,j,p,q);
                    if summ == 3
                        % if one node has 3 degree
                        if degreeofthree == 0
                            % A                            
                            if numofintra == 2
                                count_con_clique(cliquetyep,1) = count_con_clique(cliquetyep,1) + 1;
                            elseif numofintra == 1
                                count_con_clique(cliquetyep,2) = count_con_clique(cliquetyep,2) + 1;
                            elseif numofintra == 0
                                count_con_clique(cliquetyep,3) = count_con_clique(cliquetyep,3) + 1;
                            else
                                'other subtype in A'
                            end
                        else
                            % B
                            if numofintra == 1
                                count_con_clique(cliquetyep,4) = count_con_clique(cliquetyep,4) + 1;
                            else
                                'other subtype in B'
                            end
                        end
                    elseif summ == 4
                        % if one node has 3 degree
                        if degreeofthree == 0
                            % C
                            if numofintra == 2
                                count_con_clique(cliquetyep,5) = count_con_clique(cliquetyep,5) + 1;
                            elseif numofintra == 0
                                count_con_clique(cliquetyep,6) = count_con_clique(cliquetyep,6) + 1;
                            else
                                'other subtype in C'
                            end
                        else
                            % D
                            if numofintra == 2
                               count_con_clique(cliquetyep,7) = count_con_clique(cliquetyep,7) + 1;
                            elseif numofintra == 1
                                count_con_clique(cliquetyep,8) = count_con_clique(cliquetyep,8) + 1;
                            else
                                'other subtype in D'
                            end
                        end
                    elseif summ == 5
                        % E
                        if numofintra == 2
                            count_con_clique(cliquetyep,9) = count_con_clique(cliquetyep,9) + 1;
                        elseif numofintra == 1
                            count_con_clique(cliquetyep,10) = count_con_clique(cliquetyep,10) + 1;
                        else
                            'other subtype in E'
                        end
                    elseif summ == 6
                        % F
                        count_con_clique(cliquetyep,11) = count_con_clique(cliquetyep,11) + 1;
                    end
                end
            end
        end
    end
end
count_con_clique2 = count_con_clique;

% is there any node with degree of three
function bool = getdegreeofnode(structuregraph,sitesondomain1, sitesondomain2, i,j,p,q)

summ1 =  structuregraph(sitesondomain1(i), sitesondomain1(j)) ...
    + structuregraph(sitesondomain1(i), sitesondomain2(p)) ...
    + structuregraph(sitesondomain1(i), sitesondomain2(q));

summ2 =  structuregraph(sitesondomain1(j), sitesondomain1(i)) ...
    + structuregraph(sitesondomain1(j), sitesondomain2(p)) ...
    + structuregraph(sitesondomain1(j), sitesondomain2(q));

summ3 =  structuregraph(sitesondomain2(p), sitesondomain1(i)) ...
    + structuregraph(sitesondomain2(p), sitesondomain1(j)) ...
    + structuregraph(sitesondomain2(p), sitesondomain2(q));

summ4 =  structuregraph(sitesondomain2(q), sitesondomain1(i)) ...
    + structuregraph(sitesondomain2(q), sitesondomain1(j)) ...
    + structuregraph(sitesondomain2(q), sitesondomain2(p));

if summ1 == 3 || summ2 == 3 || summ3 == 3 || summ4 == 3
    bool = 1;
else
    bool = 0;
end

function numofintra = countintraedge(structuregraph, sitesondomain1, sitesondomain2, i,j,p,q)
numofintra = structuregraph(sitesondomain1(i), sitesondomain1(j)) + structuregraph(sitesondomain2(p), sitesondomain2(q));

