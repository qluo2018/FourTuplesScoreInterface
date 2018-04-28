%% This demo is a quick start for scoring the predicted interface between two chains by compare the
% local network patterns, in terms of 4 tuples, at the predicted interface with the local network
% patterns established from a data set of interfaces.

%% Thank you for using this program. Any questions about this code may send to mrqiangluo@gmail.com

%% Please cite:
%  Q Luo, R Hamer, G Reinert, CM Deane. Local Network Patterns in Protein-Protein Interfaces. PLoS One, 2012. (submitted)
%% This code is intend for acadamic use only. 

%% The current version by Qiang Luo avaliable on Feb 15, 2012.

%% Requirement:  This code runs by Matlab R2009b on Unix, and also make sure JOY is avaliable in your system
%  the subrouting CliqueCounting4forcon can be compiled to be mex
%  function to improve the speed of the counting program, see
%  CliqueCounting4forcon.m for more details

%% Get things ready for the program:
% 1. Predicted structure in PDB file;  
%    only two chains in this structure and the interface is supposed to be between these two chains
% for example, take one chain from protein A, and take the other chain from
% protein B, and name the chain to be A and B to make a new pdb file
 

%% data path and parameters setting
dataforcomplex = '';  % path for your data folder where the predicted structure is, if neccessary
temppath = ''; % path, if neccessary
resultpath = ''; % path, if neccessary
protein{1} = '1ku6';
load('stats4tuple.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare protein structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for each protein structue
% % % define the protein: name, domain, number and chain
% % % get sequence from pdb file
% % % get atom-atom distance from pdb file
% % % run Joy
% % % do statistics
% % end
%% get sequence for pdb file
'getting sequences...'
fidfasta = fopen([dataforcomplex, 'dockingcomplexallsequences.fa'],'w');
fidchain = fopen([dataforcomplex, 'dockingcomplexallchainID.txt'],'w');
fidresnum = fopen([dataforcomplex, 'dockingcomplexallresinum.txt'],'w');
for pri = 1 : size(protein,2)
    if mod(pri, 10) == 0
        pri
    end        
    clear currentprotein;
    Sequence = []; ResidueNum =[]; ChainID = [];
    currentprotein = protein{pri};
    % getpdb       
    currentprotein.structure_info = pdbread([dataforcomplex,currentprotein.names,'.pdb']);
    clear model;
    model = currentprotein.structure_info.Model;
    % get domain definition
    clear domain;
    domain.chain{1} = model.Atom(1).chainID; 
    domain.number{1}(1,1) = model.Atom(1).resSeq;
    t = 1;
    for i = 2 : size(model.Atom,2)
        key = 1;
        for j = 1 : size(domain.chain,2)
            if strcmpi(model.Atom(i).chainID, domain.chain{j})
                key = 0;
                break;
            end
        end
        if key == 1
            % new chain
            t = t + 1;
            domain.chain{t}  = model.Atom(i).chainID;
            domain.number{1}(1,2) = model.Atom(i-1).resSeq;
            domain.number{2}(1,1) = model.Atom(i).resSeq;
        end
    end
    domain.number{2}(1,2) = model.Atom(i).resSeq;    
    if t ~= 2
        protein{pri}.names
        'this protein does not exactly have two chains'
    end    
    currentprotein.aln.domain = domain;
    % get sequence   
    for i = 1 : size(model.Atom,2)
        key = 0;
        for j = 1 : 2
            for k = 1 : size(currentprotein.aln.domain.number{j},1)
                if model.Atom(i).resSeq >= currentprotein.aln.domain.number{j}(k,1) ...
                        && model.Atom(i).resSeq <= currentprotein.aln.domain.number{j}(k,2) ...
                        && model.Atom(i).chainID == currentprotein.aln.domain.chain{j}
                    key = 1;
                end
            end
        end
        if key == 1
            if size(Sequence,2) == 0
                Sequence = [Sequence, myaminolookup(model.Atom(i).resName)];
                ResidueNum = [ResidueNum, model.Atom(i).resSeq];
                ChainID = [ChainID, model.Atom(i).chainID];
            else
                if ResidueNum(size(ResidueNum,2)) == model.Atom(i).resSeq...
                        && ChainID(size(ChainID,2)) == model.Atom(i).chainID
                else
                    Sequence = [Sequence, myaminolookup(model.Atom(i).resName)];
                    ResidueNum = [ResidueNum, model.Atom(i).resSeq];
                    ChainID = [ChainID, model.Atom(i).chainID];
                end
            end
        end        
    end 
    Header = num2str(pri);
    fprintf(fidfasta, '%s ', Header);
    fprintf(fidfasta, '%s\n', Sequence);
    fprintf(fidchain, '%s\n', ChainID);
    fprintf(fidresnum, '%s\n', num2str(ResidueNum));
end
fclose(fidfasta);
fclose(fidchain);
fclose(fidresnum);

%% get atom-atom distance from pdb file
'preprocessing...'
fidfasta = fopen([dataforcomplex, 'dockingcomplexallsequences.fa'],'r');
fidchain = fopen([dataforcomplex, 'dockingcomplexallchainID.txt'],'r');
fidresnum = fopen([dataforcomplex, 'dockingcomplexallresinum.txt'],'r');
for pri = 1 : size(protein,2)    
    if mod(pri, 10) == 0
        pri
    end
    clear currentprotein;
    currentprotein = protein{pri};
    fid = fopen([temppath, currentprotein.names, '.mat']);
    if fid == -1
        % do statistics one by one
        currentprotein.gi_number = pri;
        currentprotein.alnfile_name = [temppath, 'aln',currentprotein.names,'_seqs_nr90.fasta'];
        currentprotein.infile = [temppath, currentprotein.names, '_nogap.fasta'];
        tline = fgets(fidresnum); if ~ischar(tline), break, end;
        clear str;
        str = mysplit_s(tline);
        for j = 1 : size(str,2)
            currentprotein.aln.residue.number(j) = str2num(str{j});
        end
        tline = fgets(fidchain); if ~ischar(tline), break, end;
        if isspace(tline(size(tline,2)))
            currentprotein.aln.residue.chain = tline(1:size(tline,2)-1);
        else
            currentprotein.aln.residue.chain = tline;
        end
        tline = fgets(fidfasta); if ~ischar(tline), break, end;
        clear str;
        str = mysplit(tline,' ');
        if str2num(str{1}) == pri
            currentprotein.aln.seq.Header = str{1};
            if isspace(str{2}(size(str{2},2)))
                currentprotein.aln.seq.Sequence = str{2}(1:size(str{2},2)-1);
            else
                currentprotein.aln.seq.Sequence = str{2};
            end
        else
            disp('something wrong in residue fasta file');
            pause();
        end
        fastawrite([temppath, 'aln',currentprotein.names,'_seqs_nr90.fasta'], currentprotein.aln.seq);
        clear str; clear tline;
        currentprotein.structure = pdbread([dataforcomplex,currentprotein.names,'.pdb']);
        clear model;
        if size(currentprotein.structure.Model,2) > 1
            model = currentprotein.structure.Model(1);
            currentprotein.structure.Model = [];
            currentprotein.structure.Model = model;
        else
            model = currentprotein.structure.Model;
        end
        clear domain;
        domain.chain{1} = model.Atom(1).chainID;
        domain.number{1}(1,1) = model.Atom(1).resSeq;
        t = 1;
        for i = 2 : size(model.Atom,2)
            key = 1;
            for j = 1 : size(domain.chain,2)
                if strcmpi(model.Atom(i).chainID, domain.chain{j})
                    key = 0;
                    break;
                end
            end
            if key == 1
                % new chain
                t = t + 1;
                domain.chain{t}  = model.Atom(i).chainID;
                domain.number{1}(1,2) = model.Atom(i-1).resSeq;
                domain.number{2}(1,1) = model.Atom(i).resSeq;
            end
        end
        domain.number{2}(1,2) = model.Atom(i).resSeq;
        if t ~= 2
            protein{pri}.names
            'this protein does not exactly have two chains'
        end
        currentprotein.aln.domain = domain;
        currentprotein.aln.length = size(currentprotein.aln.seq.Sequence,2);
        currentprotein.aln.locationMSA = 1 : currentprotein.aln.length;
        %             'preprocessing...'
        currentprotein = preprocess_aln_file_for_stats(currentprotein, temppath);
        save([temppath, currentprotein.names, '.mat'], 'currentprotein');
    else
        fclose(fid);
        tline = fgets(fidresnum); if ~ischar(tline), break, end;
        tline = fgets(fidchain); if ~ischar(tline), break, end;
        tline = fgets(fidfasta); if ~ischar(tline), break, end;
    end
end
fclose(fidfasta);
fclose(fidchain);
fclose(fidresnum);

%% run Joy
'joyinformation....'
skipset = []; % joy cannot run
for pri = 1 : size(protein,2)
        pri
    if ~isin(pri, skipset)
        clear currentprotein;        
        fiddata = fopen([temppath, protein{pri}.names, 'joyfeatures.mat']);
        if fiddata == -1
            load([temppath, protein{pri}.names, '.mat']);
            clear joyinformation;
            joyinformation = JoyInformation2(currentprotein, pdbpath, pdbpath2, temppath);
            EBindex = joyinformation.EBindex;
            joypsa = joyinformation.joypsa;
            EB = joyinformation.EB;
            SS = joyinformation.SS;
            secondaryindex = joyinformation.secondaryindex;
            save([temppath, protein{pri}.names, 'joyfeatures.mat'], 'EBindex','joypsa','EB','SS','secondaryindex');
            clear joyinformation;
        else
            fclose(fiddata);
        end
    end
end

% use the SA change to find the real interdomain contact from the
% candidates whose distance less than 4.5 A
% change the contact pair and contact sites in protein
% find interdomain contacts...'
asachange = [];
exceptionlist = [];
for pri = 1 : size(protein,2)    
    fid1 = fopen([temppath, protein{pri}.names, '.mat']);
    fid2 = fopen([temppath, protein{pri}.names, 'joyfeatures.mat']);
    if fid1 == -1 || fid2 == -1
        exceptionlist = [exceptionlist, pri];
    else
        fclose(fid1); fclose(fid2);
        clear currentprotein;
        load([temppath, protein{pri}.names, '.mat']);
        clear joypsa;
        load([temppath, protein{pri}.names, 'joyfeatures.mat']);
        clear result;
        result = FindInterdomainContact(currentprotein, joypsa);
        clear currentprotein;
        currentprotein = result.protein{1};
        save([temppath, currentprotein.names, '.mat'], 'currentprotein');
        asachange = [asachange,result.ASAChange];
        clear result;
    end
end

%% do stats for each protein
'stats...'
list = [1]; 
for i = 1 : size(list,2)
    pri =  list(i);
    clear currentprotein;
    load([temppath, protein{pri}.names, '.mat']);
    clear EBindex;
    load([temppath, currentprotein.names, 'joyfeatures.mat']);
    result = StatisticsWithCountNumbersForBigDatabase4nodes(currentprotein, EBindex, temppath);
    countingresult{pri} = result;
    save([resultpath, 'statsondocking.mat'], 'countingresult');
end

%% chi-square signal
% based on domain-domain interfaces
clear bgdis
clear O
clear contribution
clear E
bgdis = sum(ddi.all4clique,2) / sum(sum(ddi.all4clique,2)); 
O = sum(countingresult{1}.count_con_clique,2);
total = sum(O);
k = 210;
c = 6;
for j = 1 : 210
    E(j) = total * bgdis(j);
    contribution(j) = (O(j)-E(j))^2 / E(j);
end
contributionDDI(1,:) = contribution;
% based on homodimer interfaces
clear bgdis
clear O
clear contribution
clear E
bgdis = sum(homocomplex.all4clique,2) ./ sum(sum(homocomplex.all4clique)); 
O = sum(countingresult{1}.count_con_clique,2);
total = sum(O);
k = 210;
c = 6;
for j = 1 : 210
    E(j) = total * bgdis(j);
    contribution(j) = (O(j)-E(j))^2 / E(j);
end
contributionHomo(1,:) = contribution;
% based on heterodimer interfaces
clear bgdis
clear O
clear contribution
clear E
bgdis = sum(homocomplex.all4clique,2) ./ sum(sum(homocomplex.all4clique)); 
O = sum(countingresult{1}.count_con_clique,2);
total = sum(O);
k = 210;
c = 6;
for j = 1 : 210
    E(j) = total * bgdis(j);
    contribution(j) = (O(j)-E(j))^2 / E(j);
end
contributionHete(1,:) = contribution;
% display the signal
subplot(3,1,1)
plot(contributionDDI)
title('background on domain-domain interfaces')
subplot(3,1,2)
plot(contributionHomo)
ylabel('chi-square signal')
title('background on homodimers interfaces')
subplot(3,1,3)
plot(contributionHete)
xlabel('4-tuple type')
title('background on heterodimers interfaces')

%% score - based on DDI
reference = sum(ddi.all4clique,2) / sum(sum(ddi.all4clique,2));
for i = 1 : size(list,2)
    obs = sum(countingresult{list(i)}.count_con_clique,2);
    score1(i) = decoyscoring(obs,reference);
    names{i} = protein{list(i)}.names;
end
disp(['Score by patterns of 4-tuples at domain-domain interfaces is ', num2str(score1)]);
%% score the structures _ based on homodimer interfaces
reference = sum(homocomplex.all4clique,2) / sum(sum(homocomplex.all4clique));
for i = 1 : size(list,2)
    obs = sum(countingresult{list(i)}.count_con_clique,2);
    score3(i) = decoyscoring(obs,reference);
    names{i} = protein{list(i)}.names;
end
disp(['Score by patterns of 4-tuples at homodimer interfaces is ', num2str(score3)]);
%% score _ based on heterodimer interfaces
reference = sum(heterocomplex.all4clique,2) / sum(sum(heterocomplex.all4clique));
for i = 1 : size(list,2)
    obs = sum(countingresult{list(i)}.count_con_clique,2);
    score2(i) = decoyscoring(obs,reference);
    names{i} = protein{list(i)}.names;
end
disp(['Score by patterns of 4-tuples at heterodimer interfaces is ', num2str(score2)]);














