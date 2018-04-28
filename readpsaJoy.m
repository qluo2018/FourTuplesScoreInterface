% reading joy result .psa to have the solvent accessibility number for each
% residue
% input: .psa file, .pdf file
% output: accessibility number for each residue
function joypsa = readpsaJoy(filename, pdbfile)
clear joypsa;
fid0 = fopen(filename, 'r');
tline = [];
while 1
    tline = fgets(fid0);
    if ~ischar(tline),   break,   end
    if size(tline,2) > 6
        if strcmpi(tline(1:6), 'ACCESS')
            break;
        end
    end
end
% %     #       Res   Res   All atoms   Non P side  Polar Side  Total Side  Main Chain
% %      #       Num  type    Sum  Per.   Sum  Per.   Sum  Per.   Sum  Per.   Sum  Per.
% %      ACCESS   16   ILE     .56  1.0    .01   .0    .00   .0    .01   .0    .55  5.6

resinPDB = -10000;
t = 0;
for i = 1 : size(pdbfile.Model.Atom,2)
    if pdbfile.Model.Atom(i).resSeq ~= resinPDB  
        resinPDB = pdbfile.Model.Atom(i).resSeq;
        % split
        mystr = splittlinebyspace(tline);
        resnum = str2num(mystr{2});
        if size(mystr{3},2) == 3
            resname = mystr{3};
            signn = 0;
        else
            resname = mystr{3}(1:3);
            signn = 1;
        end
        totalside = str2num(mystr{11});
        totalsidesum = str2num(mystr{10});
        % check and put
        if resnum == resinPDB && strcmpi(resname,pdbfile.Model.Atom(i).resName)
            t = t + 1;
            joypsa(t).resnum = resnum;
            joypsa(t).resname = resname;
            joypsa(t).reschain = pdbfile.Model.Atom(i).chainID;
            joypsa(t).totalside = totalside;
            joypsa(t).atommissinginsidechain = signn;
            joypsa(t).totalsidesum = totalsidesum;
        end           
        % read again
        tline = fgets(fid0); if ~ischar(tline),   break,   end
        if size(tline,2) > 6
            if ~strcmpi(tline(1:6), 'ACCESS')
                tline = fgets(fid0); if ~ischar(tline),   break,   end
            end
        end
    end    
end
fclose(fid0);

% split the tline into words
function mystr = splittlinebyspace(tline)
% the split symbol is space
clear mystr;
mystr{1} = tline(1:6);
if ~strcmpi(mystr{1}, 'ACCESS')
    'error in reading psa'
    pause();
end
num_line = 7;
while tline(num_line) == ' '
    num_line = num_line + 1;
end
% resnum
mystr{2} = [];
while tline(num_line) ~= ' '
    mystr{2} = [mystr{2}, tline(num_line)];
    num_line = num_line + 1;
end
while tline(num_line) == ' '
    num_line = num_line + 1;
end
% resName
mystr{3} = [];
while tline(num_line) ~= ' '
    mystr{3} = [mystr{3}, tline(num_line)];
    num_line = num_line + 1;
end
while tline(num_line) == ' '
    num_line = num_line + 1;
end
% maybe a '!' 
% A `!' after the residue
%      type means that not all atoms were found in that sidechain,
%      the absolute accessibilities will be correct but the meaning
%      of the relative values is questionable, joy automatically
%      considers these residues to be accessible. 
if tline(num_line) == '!'
    mystr{3} = [mystr{3}, tline(num_line)];
    num_line = num_line + 1;
end
while tline(num_line) == ' '
    num_line = num_line + 1;
end
% numbers no more than 5 characters
for i = 4 : 13  
    mystr{i} = [];
    length = 0;
    while tline(num_line) ~= ' ' && length < 5 && num_line < size(tline,2)
        length = length + 1;
        mystr{i} = [mystr{i}, tline(num_line)];
        num_line = num_line + 1;
    end  
    while tline(num_line) == ' ' && num_line < size(tline,2)
        num_line = num_line + 1;
    end
end