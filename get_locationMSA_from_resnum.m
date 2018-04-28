%% get the corresponding locationMSA for resnum
% input: locationMSA, residue, target_res_number, target_res_chain
% input a matrix return a matrix
% if the residue is not there, location 0 will be returned
% so, need to check if the location is 0  or not before you use the result
% returned by this function

function location = get_locationMSA_from_resnum(locationMSA, residue, cpr, cpr_c)
[n,m] = size(cpr);
location1 = zeros(n,m);
for i = 1 : n
    for j = 1 : m
        for k = 1 : size(residue.number,2)
            if cpr(i,j) == residue.number(k) && cpr_c(i,j) == residue.chain(k)
                location1(i,j) = locationMSA(k);
                break;
            end
        end
    end
end
location = location1;
if size(location,1) ~= size(cpr,1)
    'error in get_locationMSA_from_resnum';
end
if size(location,2) == 0 
    location = 0;
end
        