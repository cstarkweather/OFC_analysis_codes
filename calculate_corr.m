function  [rho sig] = calculate_corr(mat1,mat2)

% num_rows = size(mat1, 1);
% corr_vals = zeros(num_rows, 1);
% 
% % Compute correlation coefficient for each row
% for i = 1:num_rows
%     r = corrcoef(mat1(i,:), mat2(i,:));
%     corr_vals(i) = r(1,2); % off-diagonal element is the correlation
% end
% 
% % To find the mean correlation
% mean_corr = mean(corr_vals);
% 
% 
% % If zero correlation falls outside these bounds, correlation is significant
% is_significant = ttest(corr_vals);
% rho = mean_corr;
% sig = is_significant;

noise_i = reshape(mat1', [], 1);
noise_j = reshape(mat2', [], 1);
[rho sig] = corr(noise_i,noise_j);

end