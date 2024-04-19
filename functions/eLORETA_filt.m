function [C, dat, filt, hasfilter, hasleadfield, hasmom, headmodel, i, keepfilter, keepleadfield, keepmom, lambda, leadfield, leadfieldopt, Nchan, Ndip, Nori, originside, origpos, rank_lf, sens, sourcemodel, varargin] = eLORETA_filt(sourcemodel, sens, headmodel, dat, C, varargin)

% Original file: ft_inverse_eloreta (Fieldtrip toolbox)
% Modified file: eLORETA_filt
% Date: April 15, 2024
% Copyright (C) 2024 Mahdi Babaei

% This file has been modified from its original version. The following 
% changes have been made:
% 1. The headmodel and the spatial filters are the same. In order to run
% the code several times and avoid the unnecessary compution of the filters
% in each iteration, calculating the eLORETA filters is seperated
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


% FT_INVERSE_ELORETA estimates the source activity using eLORETA
%
% Use as
%   [estimate] = ft_inverse_eloreta(sourcemodel, sens, headmodel, dat, cov, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'keepfilter'       = remember the spatial filter,    can be 'yes' or 'no'
%   'keepleadfield'    = remember the forward computation,  can be 'yes' or 'no'
%   'keepmom'          = remember the dipole moment,        can be 'yes' or 'no'
%   'lambda'           = scalar, regularisation parameter (default = 0.05)
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% If the dipole definition only specifies the dipole location, a rotating dipole
% (regional source) is assumed on each location. If a dipole moment is specified, its
% orientation will be used and only the strength will be fitted to the data.
%
% This implements: 
% - R.D. Pascual-Marqui; Discrete, 3D distributed, linear imaging methods of electric
%   neuronal activity. Part 1: exact, zero error localization. arXiv:0710.3341 
%   2007-October-17, http://arxiv.org/pdf/0710.3341
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Copyright (C) 2013, Marlene Boenstrup, Jan-Mathijs Schoffelen and Guido Nolte
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% if mod(nargin-5,2)
%   % the first 5 arguments are fixed, the other arguments should come in pairs
%   ft_error('invalid number of optional arguments');
% end

% get the optional input arguments, or use defaults
keepfilter      = ft_getopt(varargin, 'keepfilter', 'no');
keepmom         = ft_getopt(varargin, 'keepmom', 'yes');
keepleadfield   = ft_getopt(varargin, 'keepleadfield', 'no');
lambda          = ft_getopt(varargin, 'lambda', 0.05);

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(varargin, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(varargin, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));

% convert the yes/no arguments to the corresponding logical values
keepfilter     = istrue(keepfilter);
keepmom        = istrue(keepmom);
keepleadfield  = istrue(keepleadfield);

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(sourcemodel, 'mom');
hasleadfield  = isfield(sourcemodel, 'leadfield');
hasfilter     = isfield(sourcemodel, 'filter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(sourcemodel, 'inside')
  if hasfilter
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.filter);
  elseif hasleadfield
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.leadfield);
  else
    sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel);
  end
end

% convert to logical representation
sourcemodel = fixinside(sourcemodel);

% keep the original details on inside and outside positions
originside = sourcemodel.inside;
origpos    = sourcemodel.pos;

% select only the dipole positions inside the brain for scanning
sourcemodel.pos    = sourcemodel.pos(originside,:);
sourcemodel.inside = true(size(sourcemodel.pos,1),1);

if hasmom
  sourcemodel.mom = sourcemodel.mom(:,originside);
end

if hasfilter
  ft_info('using precomputed filters\n');
  sourcemodel.filter = sourcemodel.filter(originside);
elseif hasleadfield
%   ft_info('using precomputed leadfields\n');
  sourcemodel.leadfield = sourcemodel.leadfield(originside);
else
  ft_info('computing forward model on the fly\n');
  if hasmom
    for i=size(sourcemodel.pos,1)
      % compute the leadfield for a fixed dipole orientation
      sourcemodel.leadfield{i} = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
    end
  else
    for i=1:size(sourcemodel.pos,1)
      % compute the leadfield
      sourcemodel.leadfield{i} = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
    end
  end
end

% use existing filters, or compute them
if ~hasfilter
  % deal with reduced rank
  % check the rank of the leadfields, and project onto the lower dimensional
  % subspace if the number of columns per leadfield > rank: this to avoid
  % numerical issues in the filter computation
  rank_lf = zeros(1,size(sourcemodel.pos,1));
  for i=1:size(sourcemodel.pos,1)
    rank_lf(i) = rank(sourcemodel.leadfield{i});
  end
  if ~all(rank_lf==rank_lf(1))
    ft_error('the forward solutions have a different rank for each location, which is not supported');
  end
  if rank_lf(1)<size(sourcemodel.leadfield{1})
    ft_notice('the forward solutions have a rank of %d, but %d orientations\n',rank_lf(1),size(sourcemodel.leadfield{1},2));
    ft_notice('projecting the forward solutions on the lower dimensional subspace\n');
    for i=1:size(sourcemodel.pos,1)
      [u,s,v{i}] = svd(sourcemodel.leadfield{i}, 'econ');
      sourcemodel.leadfield{i} = sourcemodel.leadfield{i}*v{i}(:,1:rank_lf(i));
    end
  end

  % convert the leadfield into Nchan*Ndip*Nori
  [Nchan, Nori] = size(sourcemodel.leadfield{1});
  Ndip          = numel(sourcemodel.leadfield);
  leadfield     = permute(reshape(cat(2,sourcemodel.leadfield{:}),Nchan,Nori,Ndip),[1 3 2]);
    
  filt = mkfilt_eloreta_inourcode(leadfield, lambda);
  for i=1:size(sourcemodel.pos,1)
    sourcemodel.filter{i,1} = squeeze(filt(:,i,:))';
  end
end


end