function EEGdata = Select_Windows(AllData, t_start, t_finish)

% This function selects the data in a period indicated by starting and 
% ending sample. 

% INPUT:
% 1. AllData        = The EEG data of all windows
% 2. t_start        = The starting sample
% 3. t_finish       = The ending sample

% OUTPUT:
% 1. EEGdata        = The EEG data selected in a window


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



EEGdata.hdr = AllData.hdr;
EEGdata.label = AllData.label;
EEGdata.fsample = AllData.fsample;
EEGdata.trialinfo = AllData.trialinfo;
EEGdata.cfg = AllData.cfg;
for i = 1:length(AllData.time)
    EEGdata.time{1,i} = AllData.time{1,i}(t_start:t_finish);
    EEGdata.trial{1,i} = AllData.trial{1,i}(:,t_start:t_finish);
    EEGdata.sampleinfo(i,1) = AllData.sampleinfo(i,1) + t_start - 1;
    EEGdata.sampleinfo(i,2) = AllData.sampleinfo(i,1) + t_finish - 1;
    EEGdata.cfg.trl(i,1:2) = array2table([EEGdata.sampleinfo(i,1),EEGdata.sampleinfo(i,2)]);
end

end