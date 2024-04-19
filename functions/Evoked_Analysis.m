%% Source reconstruction, Evoked Data
function [dippow, pos, inside, tri, timebin, freqbin] = Evoked_Analysis(D_n, D_p, sa, winlen, stepwin, freqrng, Res_num, flag_Sp_or_Ev, flag_Channels, Selected_Chn, flag_Lambda, Lambda, pref, suf, analysis_type, index_events, eventtype, Selected_Ev, events_name, prestim, poststim)

% Applying eLORETA function on moving window
% windows. Assuming that the EEG data is already
% preprocessed.

% INPUT: 
% 1. [D_n,D_p]      = The EEG file name and its path. (Supports .edf files)
% 2. sa             = The NYC head model .mat file
% 3. winlen         = window length (seconds)
% 4. stepwin           = stepwin of moving windows forward (seconds). If stepwin is
%                     less than winlen, the windows will have overlap. If
%                     stepwin is greater than winlen, they will not have
%                     overlaps. e.g. winlen = 10 and stepwin = 1, then first
%                     window is from 0 to 8.99 and next window is from 1 to
%                     9.99. OR if winlen = 2 and stepwin = 2, first window is
%                     from 0 to 1.99 and next one is from 2 to 3.99. OR if
%                     winlen = 2 and stepwin = 3, first window is from 0 to
%                     1.99 and next window is from 3 to 4.99.
% 5. freqrng        = desired frequency range in Hz
% 6. Res_num        = desired resolution of the NYC model. 
%                     Res_num must be 2, 5, 10 or 75 which indicates
%                     how many dipoles are used in the head model.
% 7. flag_Sp_or_Ev  = shows the EEG data type
%                     which is 1 for Spont and 2 for Evoked data.
% 8. flag_Channels  = 1 for selecting a subset
%                     available channels and 0 for using all of them.
% 9. Selected_Chn   = a set containing the 
%                     number of channels e.g. 1:64 or [1,3,45,2]
%                     Each selected channel which was not available in
%                     the NYC head model, its signal will not be
%                     used in the following steps.
%                     Enter [] if the flag_Channels was set to 0.
% 10. flag_Lambda    = 1 to use a user-defined lambda
%                     and 0 to use the default lambda which is 0.05
% 11. Lambda        = a regularization term in eLoreta funtion.
%                     Enter [] if the flag_Lambda was set to 0.
% 12. [pref,suf]    = Enter all characters before and after the channel
%                     names e.g. pref = 'EEG ' and suf = '-Ref' in 'EEG Cz-Ref'
% 13. analysis_type = Shows the analysis type, 1 for frequency analysis and
%                     and 2 for timelock analysis.
% 14. index_events  = Index of events, output of Show_EEGdata_Events
%                     function.
% 15. start         = The beginning sample of events, output of Show_EEGdata_Events
% 16. Selected_Ev   = The event number you want to select. (A scalar)
% 17. events_name   = Name of the events in the data, output of Show_EEGdata_Events 
% 18. prestim       = The period (seconds) that should be selected before the trigger
% 19. poststim      = The period (seconds) that should be selected after the trigger


% OUTPUT:
% 1. dippow         = dipoles' power. M*N matrix in which M is the number
%                     of windows and N is dipoles power in each window.
% 2. pos,inside,tri = The position of dipoles, a vector indicating whether 
%                     the dipoles are in or out of the brain, the triangles           
% 3. timebin        = time bins of timelock analysis
% 4. freqbin        = freq bins of frequency analysis


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
if (flag_Sp_or_Ev  == 1)
    disp('Error: This function is for Evoked data!')
    dippow                  = [];
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading header file
cfg                         = [];
cfg.dataset                 = [D_p,D_n];
Header                      = ft_read_header(cfg.dataset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels default 
if (flag_Channels  == 0)
    Selected_Chn            = 1:Header.nChans;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lambda default
if (flag_Lambda == 0)
    Lambda                  = 0.05;
end


%############## SELECTING THE EVENT ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting the data of the selected event.
cfg.trialdef.eventtype      = eventtype;
cfg.trialdef.eventvalue     = events_name{Selected_Ev};
cfg.trialfun                = 'ft_trialfun_general';
cfg                         = ft_definetrial(cfg);
eventsidx                   = index_events{Selected_Ev};
cfg.trl(:,2)                = array2table(table2array(cfg.trl(:,2)) + poststim*Header.Fs*ones(length(eventsidx),1) - 1);
cfg.trl(:,1)                = array2table(table2array(cfg.trl(:,1)) - prestim*Header.Fs*ones(length(eventsidx),1)); 
%############## PREPARING CHANNELS ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[channel, model_channel]    = Prepare_Channels(sa, pref, suf, Selected_Chn, Header);


%############## CONVERT EEG DATA ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting the EEG file format to be used in fieldtrip. 
cfg.channel                 = channel;
cfg.continuous              = 'yes';
EEGdata                     = ft_preprocessing(cfg);


%########### PREPARING NYC HEAD MODEL ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing NYC head model with the selected resolution
Res = ['cortex',num2str(Res_num),'K'];
[dipin, grad, headmodel]    = Prepare_NYC_Head_Model(channel, model_channel, sa, Res);
pos                         = dipin.pos;
inside                      = dipin.inside;
tri                         = sa.(Res).tri;

%########### PREPARING WINDOWS ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing the beginning and ending sample of windows.
datalen                     = length(EEGdata.trial{1,1}(1,:));
windows                     = Prepare_Windows(winlen, stepwin, datalen, EEGdata.fsample);


%############## ANALYSING THE DATA ##############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AllData                     = EEGdata;
timebin                     = []; 
freqbin                     = [];
[~, ~, filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint] = eLORETA_filt(dipin, grad, headmodel, [], [], 'keepleadfield', 'yes', 'keepmom', 'no', 'lambda', Lambda);

for win = 1:length(windows(:,1))
    disp(['Window ', num2str(win), ' from ', num2str(length(windows(:,1)))])
    t_start                 = windows(win,1);
    t_finish                = windows(win,2);
    EEGdata                 = Select_Windows(AllData, t_start, t_finish);
    trlnum                  = length(EEGdata.trial);
    if (analysis_type == 1)
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
        Cf                  = prepare_freq_matrices(cfg, freq);
        freqbin(win,:)      = freq.freq;
        
        % Saving the data of separate windows
        Covf                = Cf;

        % Converting the format of the data to single window
        freq.cfg.keeptrials = 'no';
        freq.dimord         = 'chan_freq';
        freq                = rmfield(freq, {'cumsumcnt','cumtapcnt'}); 


%########### INVERSE PROBLEM SOLUTION ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse model, dipole powers computation
        dippow(win).pow     = [];
        for fr = 1:length(freqbin(win,:))
            disp(['Frequency ', num2str(fr), ' from ', num2str(length(freqbin(win,:)))])
            for i_trial = 1:trlnum
                % By commenting line 132 of Spont_Analysis function
                % and line 68 of mkfilt_eloreta function, nothing will be displayed.
                % To understand whether the code is running correcctly and 
                % the result are valid, the displayed numbers in each window should 
                % reach zero.

                Cf                              = squeeze(Covf(i_trial,:,:,fr));
                dippow_single_window            = ft_inverse_eloreta_inourcode(Cf, [], filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint);
                dippow(win).pow(i_trial,:,fr) = dippow_single_window.pow;
            end
        end
        
        
    elseif (analysis_type == 2)
    % Timelock analysis.
        timebin(win,:)          = EEGdata.time{1, 1};
        cfg                     = [];
        cfg.channel             = EEGdata.label;
        cfg.covariance          = 'yes';
        cfg.keeptrials          = 'yes';
        cfg.removemean          = 'no';
        timelock                = ft_timelockanalysis(cfg, EEGdata);
        
        % Saving the data of separate windows
        Alltrial                = timelock.trial;
        Covt                    = timelock.cov;
        
        % Converting the format of the data to single window
        timelock.cfg.keeptrials = 'no';
        timelock.dimord         = 'chan_time';
        timelock                = rmfield(timelock, {'sampleinfo','trial','trialinfo'}); 

        
%########### INVERSE PROBLEM SOLUTION ############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse model, dipole powers computation
        dippow(win).pow         = [];
%         tic
        for i_trial = 1:trlnum
            % By commenting line 132 of Spont_Analysis function
            % and line 68 of mkfilt_eloreta function, nothing will be displayed.
            % To understand whether the code is running correcctly and 
            % the result are valid, the displayed numbers in each window should 
            % reach zero.
            
            Ct                                  = squeeze(Covt(i_trial,:,:));
%             timelock.avg                        = squeeze(Alltrial(i_trial,:,:));
%             timelock.var                        = zeros(size(Alltrial(i_trial,:,:)));
%             timelock.dof                        = ones(size(Alltrial(i_trial,:,:)));
            dippow_single_window                = ft_inverse_eloreta_inourcode(Ct, [], filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint);
            dippow(win).pow(i_trial,:,:)        = dippow_single_window.pow;
        end
%         toc
    end
end 

end