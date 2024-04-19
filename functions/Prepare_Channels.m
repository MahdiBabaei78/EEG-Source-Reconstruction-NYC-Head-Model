function [channel, model_channel] = Prepare_Channels(sa, pref, suf, Selected_Chn, Header)

% The channels that are both availbale in the NYC model and the EEG dataset
% will be stored in the 'channel' variable. Usually, there are some 
% characters before and after the channels names. For example, in the 
% NYC model there is a channel name calles 'Cz'. When the EEG data is stored 
% in .edf format, this channel might have the name of 'EEG Cz-Ref'. In this 
% example, pref = 'EEG ' and suf = '-Ref'. These characters will be added to 
% the NYC channel names to make it easier to find the same channels in other 
% functions.


% INPUT: 
% 1. sa             = the NYC head model .mat file
% 2. [pref,suf]     = Enter all characters before and after the channel
%                     names e.g. pref = 'EEG ' and suf = '-Ref' in 'EEG Cz-Ref'
% 3. Selected_Chn   = a set containing the number of channels e.g. 1:64 or [1,3,45,2]
%                     Each selected channel which was not available in
%                     the NYC head model, its signal will not be
%                     used in the following steps.
%                     Enter [] if the flag_Channels was set to 0.
% 4. Header         = the EEG data header file

% OUTPUT:
% 1. channel        = The intersection of the selected channels from the
%                     EEG data and the NYC head model.
% 2. model_channel  = The NYC channels with pref and suf added to it


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







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting NYC model channels to EDF format
model_channel = sa.clab_electrodes; 
for i = 1:length(model_channel)
    model_channel{i} = [pref,model_channel{i},suf]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the user-defined channels from the data
Selected_Chn = sort(Selected_Chn);
HeaderLabel = cell(length(Selected_Chn),1);
for i = 1:length(Selected_Chn)
    HeaderLabel{i,1} = Header.label{Selected_Chn(i),1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the channels available in NYC model
channel = ft_channelselection(model_channel, HeaderLabel);
chndif = num2str(length(Selected_Chn) - length(channel));
if (chndif > 0)
disp(['Neglecting ',num2str(chndif),' channels which are not available in the NYC model.'])
for i = 1:length(HeaderLabel)
    flag_exist = 0;
    for j = 1:length(channel)
        if(strcmp(channel{j},HeaderLabel{Selected_Chn(i),1}))
            flag_exist = 1;
            break;
        end
    end
    if(flag_exist == 0)
        disp(['Channel ',num2str(i),': ',HeaderLabel{Selected_Chn(i),1}])
    end
end   
end
end
