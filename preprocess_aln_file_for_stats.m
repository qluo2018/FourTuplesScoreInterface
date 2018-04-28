% % make a new aln file for the method
% 1. put the columns for one domain together
% 2. delete the columns with too many gaps
% 3. delete the columns with too low entropy (optional)
% 4. delete the columns without structure information
% 5. find the contact sites by structure distances
% input: protein structure
% fields for protein structure: names, gi_number, alnfile_name,
% domnumfile_num, infile, structure, aln.seq, aln.domnum, aln.resnum
% output: write a file as aln for scoring system
%         seprate_site is the length of the first domain


function protein2 = preprocess_aln_file_for_stats(protein,temppath)

protein2 = protein;
%% put the columns for each domain together 
location = protein.aln.locationMSA;
residue = protein.aln.residue;
domain = protein.aln.domain;
protein2.aln.new.residue.number = [];
protein2.aln.new.residue.chain = [];
% map the domain bounds with residue number on the MSA
t = 1;
for i = 1 : size(domain.number,2)
    for j = 1 : size(domain.number{i},1)
        for k = 1 : 2
            domainmapped(t) = get_locationMSA_from_resnum(location,...
                residue, domain.number{i}(j,k), domain.chain{i});
            if domainmapped(t) == 0
                fid = fopen('domaindefinitionmodified.txt', 'a');
                fprintf(fid, '%s:  %d', protein.names, domain.number{i}(j,k));
                'error in domain re-arrangement of columns in sequence'
                'we need to modify the domain definition'
                while ~domainmapped(t) 
                    if k == 1
                        domain.number{i}(j,k) = domain.number{i}(j,k) + 1;                        
                    else
                        domain.number{i}(j,k) = domain.number{i}(j,k) - 1;
                    end
                    domainmapped(t) = get_locationMSA_from_resnum(location,...
                        residue, domain.number{i}(j,k), domain.chain{i});
                end
                fprintf(fid, '=> %d \n', domain.number{i}(j,k));
                fclose(fid);
                protein2.aln.domain = domain;                    
            end
            t = t + 1;
        end
    end
end
% rearrange the columns into domains only for two domains
length = protein.aln.length;
index_domain1 = [];
index_domain2 = [];
for i = 1 : length
    t = 1;
    for j = 1 : size(domain.number{1},1)
        if i > domainmapped(t)-1 && i < domainmapped(t+1)+1
            index_domain1 = [index_domain1, i];             
        end
        t = t + 2;
    end    
    for j = 1 : size(domain.number{2},1)
        if i > domainmapped(t)-1 && i < domainmapped(t+1)+1
            index_domain2 = [index_domain2, i];        
        end
        t = t + 2;
    end
end 
protein2.aln.new.seprate_site = size(index_domain1,2);
clear index;
index = [index_domain1, index_domain2];
protein2.aln.new.residue = get_resnum_from_locationMSA(location, residue, index);
for i = 1 : size(protein.aln.seq,1)     
    protein2.aln.new.seq(i,1).Header = protein.aln.seq(i).Header;
    protein2.aln.new.seq(i,1).Sequence = protein.aln.seq(i).Sequence(index);
end
filename = [temppath, protein.names, '.fasta'];
fastawrite(filename, protein2.aln.new.seq);
protein2.aln.new.length = size(protein2.aln.new.seq(1).Sequence,2);
protein2.aln.new.locationMSA = 1:protein2.aln.new.length;
% protein2.aln.new.residue
% pause();

%% pre-processing with the MSA: deleting the GAPs
filename_nogap = protein.infile;
protein2.aln.nogap.seq = protein2.aln.new.seq;
protein2.aln.nogap.length = protein2.aln.new.length;
protein2.aln.nogap.residue = protein2.aln.new.residue;
protein2.aln.nogap.locationMSA = 1:protein2.aln.nogap.length;
protein2.aln.nogap.seprate_site = find_seprate_site(protein2.aln.nogap.locationMSA, protein.aln.domain, protein2.aln.nogap.residue);

%% find the contact sites
threshold_contact = 4.5;
clear temp;
temp = compute_atom_distances(protein.structure, protein2.aln.new.residue);
nd = temp.distances;
nlv = temp.left_vec; % left indexes of residues
protein2.aln.nogap.distances = nd;
if size(nlv,2) ~= size(protein2.aln.nogap.residue.number,2)
    'some residue in nogap sequence doesnt have structure info.'
    pause();
end
% pause();
% contact pairs
clear cpr; clear cpr_locationMSA; clear cpr_c;
t = 1;
for i = 1 : size(nlv,2)-1
    for j = i+1 : size(nlv,2)
        if nd(i,j) < threshold_contact
            cpr(t,1:2) = [protein2.aln.new.residue.number(nlv(i)), protein2.aln.new.residue.number(nlv(j))];
            cpr_c(t,1:2) = [protein2.aln.new.residue.chain(nlv(i)), protein2.aln.new.residue.chain(nlv(j))];
            t = t + 1;
        end
    end
end
protein2.aln.contact_pairs_residue.number = cpr;
protein2.aln.contact_pairs_residue.chain = cpr_c;

% contact sites
clear csr; clear csr_c;
csr_A = []; csr_A_c = [];
csr_B = []; csr_B_c = [];
locationMSA = get_locationMSA_from_resnum(protein2.aln.new.locationMSA,...
    protein2.aln.new.residue, cpr, cpr_c);
for j = 1 : size(locationMSA,1)
    if locationMSA(i,1) == 0 || locationMSA(i,2) == 0
        disp(['error in contactSites finding in preprocess_aln_file.m', protein2.names]);
        pause();
    end
end
for i = 1 : size(cpr,1)
    if locationMSA(i,1) <= protein2.aln.new.seprate_site...
            && locationMSA(i,2) > protein2.aln.new.seprate_site
        key = 1;
        for j = 1 : size(csr_A, 2)
            if cpr(i,1) == csr_A(j) && cpr_c(i,1) == csr_A_c(j)
                key = 0;                
            end
        end
        if key == 1
            csr_A = [csr_A, cpr(i,1)];
            csr_A_c = [csr_A_c, cpr_c(i,1)];
        end    
        key = 1;
        for j = 1 : size(csr_B, 2)
            if cpr(i,2) == csr_B(j) && cpr_c(i,2) == csr_B_c(j)
                key = 0;
            end
        end
        if key == 1
            csr_B = [csr_B, cpr(i,2)];
            csr_B_c = [csr_B_c, cpr_c(i,2)];
        end
    elseif locationMSA(i,1) > protein2.aln.new.seprate_site...
            && locationMSA(i,2) <= protein2.aln.new.seprate_site
        key = 1;
        for j = 1 : size(csr_B, 2)
            if cpr(i,1) == csr_B(j)  && cpr_c(i,1) == csr_B_c(j)
                key = 0;                
            end
        end
        if key == 1
            csr_B = [csr_B, cpr(i,1)];
            csr_B_c = [csr_B_c, cpr_c(i,1)];
        end    
        key = 1;
        for j = 1 : size(csr_A, 2)
            if cpr(i,2) == csr_A(j) && cpr_c(i,2) == csr_A_c(j) 
                key = 0;
            end
        end
        if key == 1
            csr_A = [csr_A, cpr(i,2)];            
            csr_A_c = [csr_A_c, cpr_c(i,2)];
        end        
    end    
end

protein2.aln.contact_sites.Aresidue.number = csr_A;
protein2.aln.contact_sites.Aresidue.chain = csr_A_c;
protein2.aln.contact_sites.Bresidue.number = csr_B;
protein2.aln.contact_sites.Bresidue.chain = csr_B_c;
protein2.aln.contact_sites.allresidue.number = [csr_A, csr_B];
protein2.aln.contact_sites.allresidue.chain = [csr_A_c, csr_B_c];

    
