%% use the SA change to find the real interdomain contact from the
%% candidates whose distance less than 4.5 A
% change the contact pair and contact sites in protein

function result = FindInterdomainContact(protein1, Joypsa)

if ~iscell(protein1)
    protein{1} = protein1;
else
    protein = protein1;
end

clear InterContactPair; clear InterContactSite; % distance less than 4.5 and asa decreased
clear halfFalseInterContactPair; clear FalseInterContactPair; % distance less than 4.5 but asa not decreased
count_decreased = 0; count_inde = 0;

for pri = 1 : size(protein,2)
    
    InterContactPair{pri}.number = [];
    InterContactPair{pri}.chain = [];
    
    InterContactSite{pri}.Aresidue.number = [];
    InterContactSite{pri}.Aresidue.chain = [];
    InterContactSite{pri}.Bresidue.number = [];
    InterContactSite{pri}.Bresidue.chain = [];
    InterContactSite{pri}.allresidue.number = [];
    InterContactSite{pri}.allresidue.chain = [];
    
    halfFalseInterContactPair{pri}.number = [];
    halfFalseInterContactPair{pri}.chain = [];

    FalseInterContactPair{pri}.number = [];
    FalseInterContactPair{pri}.chain = [];
    
    
    clear locationMSA; clear residue;
    locationMSA = protein{pri}.aln.new.locationMSA;
    residue = protein{pri}.aln.new.residue;
    seprate_site = protein{pri}.aln.new.seprate_site;
    % go through the contact pairs
    location = get_locationMSA_from_resnum(locationMSA, residue, ...
        protein{pri}.aln.contact_pairs_residue.number,...
        protein{pri}.aln.contact_pairs_residue.chain);    
    for i = 1 : size(location,1)
        if location(i,1) ~= 0 && location(i,2) ~= 0      
            % find inter-domain contact pairs
            if location(i,1) <= seprate_site && location(i,2) > seprate_site ...
                    || location(i,2) <= seprate_site && location(i,1) > seprate_site 
                % find if solvent accessibilities of these two sites are
                % all decreased
                % and they must both be on the surface
                psachange1 = findPSAbyResNum(protein{pri}.aln.contact_pairs_residue.number(i,1), ...
                    protein{pri}.aln.contact_pairs_residue.chain(i,1), Joypsa{pri});
                psachange2 = findPSAbyResNum(protein{pri}.aln.contact_pairs_residue.number(i,2), ...
                    protein{pri}.aln.contact_pairs_residue.chain(i,2), Joypsa{pri});
                psa1 = findPSATotalSidebyResNum(protein{pri}.aln.contact_pairs_residue.number(i,1), ...
                    protein{pri}.aln.contact_pairs_residue.chain(i,1), Joypsa{pri});
                psa2 = findPSATotalSidebyResNum(protein{pri}.aln.contact_pairs_residue.number(i,2), ...
                    protein{pri}.aln.contact_pairs_residue.chain(i,2), Joypsa{pri});               
                if  psachange1 < 0  &&  psachange2 < 0 && psa1 > 7 && psa2 > 7
                    InterContactPair{pri}.number = [InterContactPair{pri}.number; protein{pri}.aln.contact_pairs_residue.number(i,:)];
                    InterContactPair{pri}.chain = [InterContactPair{pri}.chain; protein{pri}.aln.contact_pairs_residue.chain(i,:)];
                    count_decreased =  count_decreased + 1;
                else
                    halfFalseInterContactPair{pri}.number = [halfFalseInterContactPair{pri}.number; protein{pri}.aln.contact_pairs_residue.number(i,:)];
                    halfFalseInterContactPair{pri}.chain = [halfFalseInterContactPair{pri}.chain; protein{pri}.aln.contact_pairs_residue.chain(i,:)];
                    if psa1 > 7 && psa2 > 7
                        if psachange1 < 0 || psachange2 < 0 
                            count_inde = count_inde + 1;
                        end
                    end
                end
            end
        else
            % if there is 0 in location, there must be something wrong
            'wrong in residue identification in FindInterdomainContact'
        end
    end
        
    %% change the contact pairs and sites in protein structure 
    % pair
    protein{pri}.aln.contact_pairs_residue.number = []; protein{pri}.aln.contact_pairs_residue.chain = [];
    protein{pri}.aln.contact_pairs_residue.number = InterContactPair{pri}.number;
    protein{pri}.aln.contact_pairs_residue.chain = InterContactPair{pri}.chain;
    
    % site
    % contact sites
    for i =  1 : size(InterContactPair{pri}.number, 1) 
        for p = 1 : 2
            keyy = 0;
            % A
            for k = 1 : size(protein{pri}.aln.domain.number{1},1)
                if InterContactPair{pri}.number(i,p) >= protein{pri}.aln.domain.number{1}(k,1) ...
                        && InterContactPair{pri}.number(i,p) <=  protein{pri}.aln.domain.number{1}(k,2) ...
                        && strcmpi(InterContactPair{pri}.chain(i,p), protein{pri}.aln.domain.chain{1})     
                    keyy = 1;
                    for j = 1 : size(InterContactSite{pri}.Aresidue.number, 2)
                        % if this site is not currently in
                        if InterContactSite{pri}.Aresidue.number(j) == InterContactPair{pri}.number(i,p)
                            keyy = 0;
                            break;
                        end
                    end
                end
            end
            if keyy == 0
                % B
                for k = 1 : size(protein{pri}.aln.domain.number{2},1)
                    if InterContactPair{pri}.number(i,p) >= protein{pri}.aln.domain.number{2}(k,1) ...
                            && InterContactPair{pri}.number(i,p) <=  protein{pri}.aln.domain.number{2}(k,2) ...
                            && strcmpi(InterContactPair{pri}.chain(i,p), protein{pri}.aln.domain.chain{2})
                        keyy = 2;
                        for j = 1 : size(InterContactSite{pri}.Bresidue.number, 2)
                            % if this site is not currently in
                            if InterContactSite{pri}.Bresidue.number(j) == InterContactPair{pri}.number(i,p)
                                keyy = 0;
                                break;
                            end
                        end
                    end
                end
            end  
            if keyy == 1
                InterContactSite{pri}.Aresidue.number = ...
                    [InterContactSite{pri}.Aresidue.number, InterContactPair{pri}.number(i,p)];
                InterContactSite{pri}.Aresidue.chain = ...
                    [InterContactSite{pri}.Aresidue.chain, InterContactPair{pri}.chain(i,p)];
            elseif keyy == 2
                InterContactSite{pri}.Bresidue.number = ...
                    [InterContactSite{pri}.Bresidue.number, InterContactPair{pri}.number(i,p)];
                InterContactSite{pri}.Bresidue.chain = ...
                    [InterContactSite{pri}.Bresidue.chain, InterContactPair{pri}.chain(i,p)];
            else
                %'dont do anyting'
            end 
        end
    end
    % 
    InterContactSite{pri}.allresidue.number = [InterContactSite{pri}.allresidue.number,...
        InterContactSite{pri}.Aresidue.number, InterContactSite{pri}.Bresidue.number];
    InterContactSite{pri}.allresidue.chain = [InterContactSite{pri}.allresidue.chain,...
        InterContactSite{pri}.Aresidue.chain, InterContactSite{pri}.Bresidue.chain];
    % site
    protein{pri}.aln.contact_sites.Aresidue.number = []; protein{pri}.aln.contact_sites.Aresidue.chain = []; 
    protein{pri}.aln.contact_sites.Bresidue.number = []; protein{pri}.aln.contact_sites.Bresidue.chain = [];
    protein{pri}.aln.contact_sites.allresidue.number = []; protein{pri}.aln.contact_sites.allresidue.chain = []; 
    protein{pri}.aln.contact_sites.Aresidue.number = InterContactSite{pri}.Aresidue.number;    
    protein{pri}.aln.contact_sites.Aresidue.chain = InterContactSite{pri}.Aresidue.chain;
    protein{pri}.aln.contact_sites.Bresidue.number = InterContactSite{pri}.Bresidue.number;    
    protein{pri}.aln.contact_sites.Bresidue.chain = InterContactSite{pri}.Bresidue.chain;
    protein{pri}.aln.contact_sites.allresidue.number = InterContactSite{pri}.allresidue.number;
    protein{pri}.aln.contact_sites.allresidue.chain = InterContactSite{pri}.allresidue.chain;
end

%% also find the total change of ASA on the interface
for pri = 1 : size(protein,2)
    ASAChange(pri) = 0;
    currentresidue.number = protein{pri}.aln.contact_sites.allresidue.number;
    currentresidue.chain = protein{pri}.aln.contact_sites.allresidue.chain;
    for j = 1 : size(currentresidue.number,2)
        psachange = findjoypsaforresidue(Joypsa{pri},j,currentresidue);
        ASAChange(pri) = ASAChange(pri) + psachange.sachangesum;
    end    
end
result.halfFalseInterContactPair = halfFalseInterContactPair;
result.ASAChange = ASAChange;
result.protein = protein;
x = ['SA decreased pairs: ', num2str(count_decreased),'SA decreased and increased: ', num2str(count_inde), 'In Total ', num2str(count_decreased+count_inde)]




