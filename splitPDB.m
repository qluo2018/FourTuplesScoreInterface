%% split the pdb file according to domain informations
% input: protein structure
% output: splited pdb files for domains

function splitPDB(datapath, protein)

domain = protein.aln.domain;
for i = 1 : size(domain.number,2)
    fid_domain{i} =  fopen([datapath, protein.names,'Domain', num2str(i), '.pdb'], 'w');    
end

fid0 = fopen([datapath, protein.names,'.pdb'], 'r');
tline = [];
while 1
    tline = fgets(fid0);
    if ~ischar(tline),   break,   end
    if size(tline,2) > 4
        if strcmpi(tline(1:4), 'ATOM')
            break;
        end
    end
end
for i = 1 : size(protein.structure.Model.Atom,2)
    fid = 0;
    for j = 1 : size(domain.number,2)
        % go through domains
        for k = 1 : size(domain.number{j},1)
            % go through regions of each domain
            if protein.structure.Model.Atom(i).resSeq >= domain.number{j}(k,1) ...
                    && protein.structure.Model.Atom(i).resSeq <= domain.number{j}(k,2) ...
                    && strcmpi(protein.structure.Model.Atom(i).chainID, domain.chain{j})
                % read and write
                fid = fid_domain{j};
            end
        end            
    end    
    if fid ~= 0
        % read
        fprintf(fid,'%s',tline);
%         pause();
        tline = fgets(fid0); if ~ischar(tline),   break,   end
        if size(tline,2) > 4
            if ~strcmpi(tline(1:4), 'ATOM')
                tline = fgets(fid0); if ~ischar(tline),   break,   end
            end
        end
    else
        tline = fgets(fid0); if ~ischar(tline),   break,   end
        if size(tline,2) > 4
            if ~strcmpi(tline(1:4), 'ATOM')
                tline = fgets(fid0); if ~ischar(tline),   break,   end
            end
        end
    end
end
fclose(fid0);
for i = 1 : size(domain.number,2)
    fclose(fid_domain{i});   
end


        
        
        
        
        