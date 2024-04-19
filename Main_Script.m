% This script is a conventional way to use the provided
% functions in the EEG-Source-Reconstruction repository.
% Run each section in order and seperately.

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

%%  Input Phase, Files and toolboxes
clc
clear
initpath = uigetdir([],'Specify which folder you want the folder searches to start from');
[D_n,D_p] = uigetfile(initpath,'Select the EEG data file (.eeg, .edf)');
functions_path = uigetdir(initpath,'Select the functions folder');
[Model_name,Model_path] = uigetfile(initpath,'Select the NYC head model (sa_nyhead.mat)');
Fieldtrip_path = uigetdir(initpath,'Select the Fieldtrip folder');
% eeglab_path = uigetdir(initpath,'Select the eeglab folder');
cd(Fieldtrip_path)
addpath('trialfun','utilities','template','statfun','src','specest','realtime','qsub','preproc','plotting','inverse','forward','fileio','contrib','connectivity','bin','privatee')
addpath(Fieldtrip_path)
% cd(eeglab_path)
% addpath('plugins\Biosig3.8.3\biosig\t200_FileAccess','plugins\Biosig3.8.3\biosig\doc','plugins\Biosig3.8.3\biosig\t250_ArtifactPreProcessingQualityControl',genpath('functions'))
% addpath(eeglab_path)
addpath(functions_path)
cd(initpath)
addpath(Model_path)
load([Model_path,Model_name]);

%% Spont Data, Freq Analysis, non-overlapping win
Show_EEGdata_Channels(D_n, D_p);
%%
winlen = 3;
overlap = 0;
freqrng = [8,12];
Res_num = 5;
flag_Sp_or_Ev = 1;
flag_Channels = 1;
Selected_Chn = 1:62;
flag_Lambda = 0;
Lambda = [];
pref = [];
suf = [];
[dippow, pos, inside, tri, freqbin] = Spont_Analysis(D_n, D_p, sa, winlen, overlap, freqrng, Res_num, flag_Sp_or_Ev, flag_Channels, Selected_Chn, flag_Lambda, Lambda, pref, suf);
%% Evoked Data, Freq or Timelocked Analysis 
Show_EEGdata_Channels(D_n, D_p);
prestim             = 1;
poststim            = 1;
[eventtype, events_name, index_events] = Show_EEGdata_Events('G:\Daneshgah\Internship Juel\Last Code\EEG files','SUB_1001_noload_p3b.edf', prestim, poststim);
%%
winlen              = 1;
stepwin             = 0.05;
freqrng             = [1,45];
Res_num             = 1;
flag_Sp_or_Ev       = 2;
flag_Channels       = 1;
Selected_Chn        = 1:62;
flag_Lambda         = 0;
Lambda              = [];
pref                = 'EEG ';
suf                 = '-Ref';
analysis_type       = 2;
Selected_Ev         = 2;
[dippow, pos, inside, tri, timebin, freqbin] = Evoked_Analysis(D_n, D_p, sa, winlen, stepwin, freqrng, Res_num, flag_Sp_or_Ev, flag_Channels, Selected_Chn, flag_Lambda, Lambda, pref, suf, analysis_type, index_events, eventtype, Selected_Ev, events_name, prestim, poststim);
















