%% Generate ordered vs shuffled choice sequences
Nblocks = 100;
offset1 = 1;
offset2 = 4;
blocklen = 20; %max len of block

blocks1 = zeros(Nblocks, blocklen);
blocks2 = zeros(Nblocks, blocklen);
blocks1(:, 1:offset1) = 1;
blocks2(:, 1:offset2) = 1;

blockarr = [blocks1; blocks2];
blockarr_ori = blockarr;

% Now shuffle within each column
for i = 1:blocklen
    col = blockarr(:,i);
    order = randperm(Nblocks * 2);
    colshuffled = col(order);
    blockarr(:,i) = colshuffled;
    
end


%% Perform regression



