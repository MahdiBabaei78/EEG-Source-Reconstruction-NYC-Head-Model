function [eventtype, events_name, index_events] = Show_EEGdata_Events(D_n, D_p, prestim, poststim)

% Show_EEGdata_Events reads data from a EDF file with channels that have a different
% sampling rates. It extracts the events and shows a list of them in the
% command window. 
%
% Use as
%   [eventtype, events_name, index_events] = Show_EEGdata_Events(D_n, D_p, prestim, poststim)
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






[EEG,EEG.event] = edf2fieldtrip('G:\Daneshgah\Internship Juel\Last Code\EEG files\SUB_1001_noload_p3b.edf');


event = EEG.event;
if isempty(event)
    disp('There is no event available in the data!')
    return;
end
eventtype = event(1,1).type;


for i = 1:length(EEG.event)
    all_events{i} = EEG.event(i).value;
    start(i) = EEG.event(i).sample;
%     finish(i) = EEG.event(i).latency + EEG.event(i).duration*EEG.srate - 1;
%     duration(i) = EEG.event(i).duration;
end

events_name{1} = all_events{1};
counter = 1;
for i = 2:length(EEG.event)
    ev = EEG.event(i).value;
    flag_exist = 0;
    for j = 1:counter
        if (strcmp(ev,events_name{j}))
           flag_exist = 1;
           break;
        end
    end
    if (flag_exist == 0)
        counter = counter + 1;
        events_name{counter} = ev;
    end
end

disp('Available Events: ')
for j = 1:length(events_name)
    
    disp(['Event ',num2str(j),': ',events_name{j}])
    index_events{j} = [];
    for i = 1:length(EEG.event) 
        if (strcmp(events_name{j},EEG.event(i).type))
           index_events{j} = [index_events{j},i]; 
        end
    end
end
disp(['Trial length is ', num2str(prestim+poststim),' seconds. Window length must be less than this value!'])

end