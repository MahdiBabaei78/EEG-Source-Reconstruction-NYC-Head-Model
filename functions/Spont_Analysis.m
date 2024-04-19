%% Source reconstruction, Spont multiple windows
function [dippow, pos, inside, tri, freqbin] = Spont_Analysis(D_n, D_p, sa, winlen, overlap, freqrng, Res_num, flag_Sp_or_Ev, flag_Channels, Selected_Chn, flag_Lambda, Lambda, pref, suf)

% Applying eLoreta function on multiple
% non-overlapping windows. Assuming that the EEG
% data is already preprocessed.

% INPUT: 
% 1. [D_n,D_p]      = The EEG file name and its path. (Supports .edf files)
% 2. sa             = The NYC head model .mat file
% 3. winlen         = window length in seconds
% 4. overlap        = number between 0 and 1 (exclusive) specifying the 
%                     fraction of overlap between snippets (0 = no overlap)
% 4. freqrng        = desired frequency range in Hz
% 5. Res_num        = desired resolution of the NYC model. 
%                     Res_num must be 2, 5, 10 or 75 which indicates
%                     how many dipoles are used in the head model.
% 6. flag_Sp_or_Ev  = shows the EEG data type
%                     which is 1 for Spont and 2 for Evoked data.
% 7. flag_Channels  = 1 for selecting a subset
%                     available channels and 0 for using all of them.
% 8. Selected_Chn   = a set containing the 
%                     number of channels e.g. 1:64 or [1,3,45,2]
%                     Each selected channel which was not available in
%                     the NYC head model, its signal will not be
%                     used in the following steps.
%                     Enter [] if the flag_Channels was set to 0.
% 9. flag_Lambda    = 1 to use a user-defined lambda
%                     and 0 to use the default lambda which is 0.05
% 10. Lambda        = a regularization term in eLoreta funtion.
%                     Enter [] if the flag_Lambda was set to 0.
% 11. pref and suf  = Enter all characters before and after the channel
%                     names e.g. pref = 'EEG ' and suf = '-Ref' in 'EEG Cz-Ref'

% OUTPUT:
% 1. dippow         = dipoles' power. M*N matrix in which M is the number
%                     of windows and N is dipoles power in each window.
% 2. pos,inside,tri = The position of dipoles, a vector indicating whether 
%                     the dipoles are in or out of the brain, the triangles           
% 3. freqbin        = freq bins of frequency analysis
% 

% Copyright (C) 2024-, Mahdi Babaei and Bjorn Erik Juel
%
% This file is part of EEG Source Reconstruction - NYC head model, see 
% https://github.com/MahdiBabaei78/EEG-Source-Reconstruction-NYC-Head-Model
% for the documentation and details.
%
%    EEG-Source-Reconstruction-NYC-Head-Model is a free software: 
%    you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EEG-Source-Reconstruction-NYC-Head-Model is distributed in the hope
%    that it will be useful, but WITHOUT ANY WARRANTY; without even the
%    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
%    PURPOSE.  See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EEG-Source-Reconstruction-NYC-Head-Model. If not, 
%    see <http://www.gnu.org/licenses/>.
%
% $Id$



%######### CHECKING ERRORS AND DEFAULTS #########
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function supports Spont data analysis.
if (flag_Sp_or_Ev  == 2)
    disp('Error: This function is for spontaneous data!')
    dippow = [];
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading header file
cfg = [];
cfg.dataset = [D_p,D_n];
Header = ft_read_header(cfg.dataset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels default 
if (flag_Channels  == 0)
    Selected_Chn = 1:Header.nChans;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lambda default
if (flag_Lambda == 0)
    Lambda = 0.05;
end


%############## TRIAL DEFINITION ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining non-overlapping windows of the signal as trials.
cfg.trialdef.triallength    = winlen;
cfg.trialdef.overlap        = overlap;
cfg.trialfun                = 'ft_trialfun_general';
cfg                         = ft_definetrial(cfg);
winnum                      = length(cfg.trl(:,1));

%############## PREPARING CHANNELS ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[channel, model_channel] = Prepare_Channels(sa, pref, suf, Selected_Chn, Header);


%############## CONVERT EEG DATA ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting the EEG file format to be used in fieldtrip. 
cfg.channel     = channel;
cfg.continuous  = 'yes';
EEGdata         = ft_preprocessing(cfg);


%########### PREPARING NYC HEAD MODEL ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing NYC head model with the selected resolution
Res = ['cortex',num2str(Res_num),'K'];
[dipin, grad, headmodel]    = Prepare_NYC_Head_Model(channel, model_channel, sa, Res);
pos                         = dipin.pos;
inside                      = dipin.inside;
tri                         = sa.(Res).tri;


%############## ANALYSING THE DATA ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency analysis with multitaper.
cfg                 = [];
cfg.method          = 'mtmfft';
cfg.output          = 'powandcsd';
cfg.channel         = EEGdata.label;
cfg.taper           = 'hanning';
cfg.foilim          = freqrng;
cfg.pad             = 'nextpow2';
cfg.tapsmofrq       = 1;
cfg.keeptrials      = 'yes';
freq                = ft_freqanalysis(cfg, EEGdata);
Cf    = prepare_freq_matrices(cfg, freq);

% Saving the data of separate windows
Covf                = Cf;

% Converting the format of the data to signle window
freq.cfg.keeptrials = 'no';
freq.dimord         = 'chan_freq';
freq                = rmfield(freq, {'cumsumcnt','cumtapcnt'}); 


%########### INVERSE PROBLEM SOLUTION ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse model, dipole powers computation
dippow.inside                   = [];
dippow.pos                      = [];
dippow.pow                      = [];
dippow.ori                      = [];
dippow.leadfield                = [];
[~, ~, filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint] = eLORETA_filt(dipin, grad, headmodel, [], [], 'keepleadfield', 'yes', 'keepmom', 'no', 'lambda', Lambda);
freqbin                         = freq.freq;
for fr = 1:length(freqbin)
    disp(['Frequency ', num2str(fr), ' from ', num2str(length(freqbin))])
    for i_win = 1:winnum
        % By commenting line 132 of Spont_Analysis function
        % and line 68 of mkfilt_eloreta function, nothing will be displayed.
        % To understand whether the code is running correcctly and 
        % the result are valid, the displayed numbers in each window should 
        % reach zero.
        disp(['Window ', num2str(i_win), ' from ', num2str(winnum)])
        Cf                       = squeeze(Covf(i_win,:,:,fr));
        dippow_single_window     = ft_inverse_eloreta_inourcode(Cf, [], filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint);
        dippow.pow(i_win,:,fr) = dippow_single_window.pow;
    end
end