function score = decoyscoring(obs,reference)
% qpo = quantile(reference, 0.975); % positive
% qne = quantile(reference, 0.025); % negative
% score = 0;
% for i = 1 : size(reference,1)
%     if reference(i) < qne
%         score = score - obs(i);
%     elseif reference(i) > qpo
%         score = score + obs(i);        
%     end        
% end

total = sum(obs);
score = 0;
for j = 1 : size(reference,1)
    E(j) = total * reference(j);
    if E(j) ~= 0
        score = score + (obs(j)-E(j))^2 / E(j);
    else
        score = score + (obs(j)-E(j))^2 ;
    end    
end
score = -score;
