function [w] = TRCA(eeg)
% -------------------------------------------------------------------
% Task-related component analysis[1].extract task-related components 
% by maximizing the reproducibility during the task period.
%
% Input:         
%   eeg : input eeg training data
%         (# channels, # points, # trials)
% Output:  
%   W: w means eigenvectors corresponding to max eigenvalues
%
% Reference:
%   [1] H. Tanaka, T. Katura, H. Sato,
%       "Task-related component analysis for functional neuroimaging and 
%        application to near-infrared spectroscopy data",
%       NeuroImage, vol. 64, pp. 308-327, 2013.
% --------------------------------------------------------------------
[num_chans, num_smpls, num_trials]  = size(eeg);
S = zeros(num_chans);
for trial_i = 1:1:num_trials-1
    x1 = squeeze(eeg(:,:,trial_i));
    x1 = bsxfun(@minus, x1, mean(x1,2));
    for trial_j = trial_i+1:1:num_trials
        x2 = squeeze(eeg(:,:,trial_j));
        x2 = bsxfun(@minus, x2, mean(x2,2));
        S = S + x1*x2' + x2*x1';
    end % end trial_j
end % end trial_i
UX = reshape(eeg, num_chans, num_smpls*num_trials);
UX = bsxfun(@minus, UX, mean(UX,2));
Q = UX*UX';
[V,~] = eigs(S, Q);
w = V(:,1);
end %end function
