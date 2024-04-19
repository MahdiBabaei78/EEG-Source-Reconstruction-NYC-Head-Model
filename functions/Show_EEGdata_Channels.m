function Show_EEGdata_Channels(D_n, D_p)

% Show_EEGdata_Channels reads data from an EDF file with channels.
% It extracts the channels and shows a list of them in the
% command window. 
%
% Use as
%    Show_EEGdata_Channels(D_n, D_p)
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


cfg = [];
cfg.dataset = [D_p,D_n];
Header = ft_read_header(cfg.dataset);
disp('Here are the available channels in your data: ')
for i = 1:Header.nChans
    disp(['Channel ',num2str(i),': ',Header.label{i,1}])
end

end