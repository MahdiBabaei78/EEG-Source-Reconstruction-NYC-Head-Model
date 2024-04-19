function [estimate] = ft_inverse_eloreta_inourcode(C, dat, filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargint)

% Original file: ft_inverse_eloreta (Fieldtrip toolbox)
% Modified file: ft_inverse_eloreta_inourcode
% Date: April 15, 2024
% Copyright (C) 2024 Mahdi Babaei

% This file has been modified from its original version. The following 
% changes have been made:
% 1. The headmodel and the spatial filters are the same. In order to run
% the code several times and avoid the unnecessary compution of the filters
% in each iteration, calculating the eLORETA objective function is seperated
% from the rest of the code.

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




% get the power
siz_C  = size(C); % C can have both a freq and time dimension
sourcemodel.pow = zeros([size(sourcemodel.pos,1),1]);
for i=1:size(sourcemodel.pos,1)
    csd               = sourcemodel.filter{i}*C(:,:)*sourcemodel.filter{i}';
    [u,s,v]           = svd(real(csd));
    sourcemodel.pow(i)    = s(1);
end
% reassign the estimated values over the inside and outside grid positions
estimate.inside  = originside;
estimate.pos     = origpos;
estimate.pow( originside,:,:) = sourcemodel.pow;
estimate.pow(~originside,:,:) = nan;
end