function windows = Prepare_Windows(winlen, stepwin, datalen, fs);

% Prepare the windows to select from the data

% INPUT: 
% 1. winlen           = window length in seconds
% 2. stepwin          = stepwin of moving windows forward (seconds). If stepwin is
%                     less than winlen, the windows will have overlap. If
%                     stepwin is greater than winlen, they will not have
%                     overlaps. e.g. winlen = 10 and stepwin = 1, then first
%                     window is from 0 to 8.99 and next window is from 1 to
%                     9.99. OR if winlen = 2 and stepwin = 2, first window is
%                     from 0 to 1.99 and next one is from 2 to 3.99. OR if
%                     winlen = 2 and stepwin = 3, first window is from 0 to
%                     1.99 and next window is from 3 to 4.99.
% 3. datalen          = The length of the trials of the original EEG data.
% 4. fs               = Sampling frequency

% OUTPUT:
% 1. windows          = a N*2 matrix, indicating the beginnig and ending
%                       sample of N windows


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





stepwin = stepwin*fs;
winlen = winlen*fs;
windows = transpose((0:floor((datalen-winlen)/stepwin)+1)*stepwin)+1;
windows = [windows, windows+winlen-1];
if (mod(datalen-winlen,stepwin) ~= 0)
    windows(end,2) = datalen;
    windows(end,1) = datalen-winlen+1;
else
    windows(end,:) = [];
end
end