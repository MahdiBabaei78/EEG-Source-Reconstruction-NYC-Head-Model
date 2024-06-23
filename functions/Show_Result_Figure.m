function Show_Result_Figure(dippow, pos, inside, tri)

% This function shows the results of the spontaneous data source 
% reconstruction in one figure. It cannot be used to 
% plot the figure of the ouput of the Spont_Analysis function directly.
% For example, an average of the dipoles' powers can be calculated 
% across all trials and frequency bins. Then, the resulting vector can 
% be plotted.

% INPUT: 
% 1. dippow         = a "vector" containing the power of each dipole (e.g. 5004*1 matrix)
% 2. pos,inside,tri = The position of dipoles, a vector indicating whether 
%                     the dipoles are in or out of the brain, the triangles 

% OUTPUT:
% A figure showing the results on the NYC brain




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


temp = dippow;
dippow = [];
dippow.pow = temp;
max_pow            = max(dippow.pow);
min_pow            = min(dippow.pow);
dippow.inside      = inside;
dippow.pos         = pos;
dippow.tri         = tri;
cfg                = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
max_op             = max_pow/2;
min_op             = min_pow/2;
cfg.funcolorlim    = [min_pow max_pow];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [min_op max_op];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'Surface_NYC.mat'; 
cfg.surfdownsample = 10;                % downsample to speed up processing
ft_sourceplot(cfg, dippow);

view ([-157 11])                        % rotate the object in the view

end
