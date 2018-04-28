%% call joy assign secondary structures and exposed/buried properties
function result = JoyInformation2(protein1, datapath, datapath2, resultpath)

if ~iscell(protein1)
    protein{1} = protein1;
else
    protein = protein1;
end

joycannot = [];
for i = 1 : size(protein,2)
     % call joy for each protein
     [s,w] = runUnix(datapath, ['psa ', datapath, protein{i}.names, '.pdb']);
     % split the pdb file for different domains
     splitPDB(datapath2, protein{i});
     EBindex{i} = zeros(1, protein{i}.aln.nogap.length);
     EB{i}.residue.num = [];
     EB{i}.residue.chain = [];
     EB{i}.character = [];
     for j = 1 : size(protein{i}.aln.domain.number,2)
          % run Joy for each domain
         [s,w] = runUnix(datapath, ['psa ', datapath, protein{i}.names, 'Domain', num2str(j), '.pdb']);
     end 
     % psa for all protein
     joypsa{i} = readpsaJoy([datapath2, protein{i}.names, '.psa'], protein{i}.structure);
     % psa for domain
     for j = 1 : size(protein{i}.aln.domain.number,2)
         clear structure_info;
         structure_info = pdbread([datapath2, protein{i}.names, 'Domain', num2str(j),'.pdb']);
         clear domainjoypsa;         
         domainjoypsa = readpsaJoy([datapath2, protein{i}.names, 'Domain', num2str(j), '.psa'], structure_info);         
         for kd = 1 : size(domainjoypsa,2)
             for ka = 1 : size(joypsa{i},2)
                 if domainjoypsa(kd).resnum == joypsa{i}(ka).resnum ...
                         && strcmpi(domainjoypsa(kd).resname, joypsa{i}(ka).resname) ...
                         && strcmpi(domainjoypsa(kd).reschain, joypsa{i}(ka).reschain)
                     joypsa{i}(ka).splitedtotalside = domainjoypsa(kd).totalside;
                     joypsa{i}(ka).splitedtotalsidesum = domainjoypsa(kd).totalsidesum;
                     joypsa{i}(ka).sachange = joypsa{i}(ka).totalside - joypsa{i}(ka).splitedtotalside;
                     joypsa{i}(ka).sachangesum = joypsa{i}(ka).totalsidesum - joypsa{i}(ka).splitedtotalsidesum;
                     break;
                 end
             end
         end
     end   
end

%% consistency check between EBindex and joypsa
for i = 1 : size(protein,2)
    clear residue; clear locationMSA;
    residue = protein{i}.aln.nogap.residue;
    locationMSA = protein{i}.aln.nogap.locationMSA;
    EBindex{i} = ones(1,protein{i}.aln.nogap.length);
    for j = 1 : size(residue.number,2)
        for k = 1 : size(joypsa{i},2)
            if residue.number(j) == joypsa{i}(k).resnum ...
                    && strcmpi(residue.chain(j), joypsa{i}(k).reschain)
                if joypsa{i}(k).splitedtotalside > 7
                    EBindex{i}(locationMSA(j)) = 2;                     
                    break;
                end                 
            end               
        end
    end   
end

result.EB = EB; result.SS = []; 
result.EBindex = EBindex; 
result.secondaryindex = [];
result.joypsa = joypsa;
