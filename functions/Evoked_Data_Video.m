% Original file: ft_sourceplot (Fieldtrip toolbox)
% Modified file: Evoked_Data_Video
% Date: April 15, 2024
% Copyright (C) 2024 Mahdi Babaei

% To visualize the results on the NYC brain using the uifigure and
% update the figure when the time slider is moved, several changes must
% have been made.
% This file has been modified from its original version. Overall, the
% modifications are as follows:

% 1. The figure handle should be given as an input to other functions
% because the gcf function will not work to obtain the UI figure handle.
% 2. The visualization parameters are specified
% 3. The UI figure, UI axes and the slider is created.
% 4. The figure updating function is implemented to work interactively and 
%    updates the plot when the slider is moved.
% 5. The script is implemented somehow to plot the brain surface once and 
%    to update the color of the dipoles as the slider moves. This will make
%    the figure to be smooth and fast because it skips unnecessary plotting.   


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



%% Input phase 
winlen = 0.006;
stepwin = 0.006;
datalen = 2000;
fsample = 1000;
windows = Prepare_Windows(winlen, stepwin, datalen, fsample);
timebin = windows*1000/fsample;
reldippow = reldippow_stim49_win6_removedmean;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear trialnum assign Ce9dei2ZOo_argin Ce9dei2ZOo_debug Ce9dei2ZOo_ds Ce9dei2ZOo_funname Ce9dei2ZOo_is Ce9dei2ZOo_ns Ce9dei2ZOo_ws cfg dim dimord dimtok doimage fcolmax fcolmin FjmoT6aA_current_ft_toplevel FjmoT6aA_ft_path FjmoT6aA_ft_ver FjmoT6aA_highest_ft FjmoT6aA_k FjmoT6aA_lowest_ft FjmoT6aA_stack ft_default ft_revision ftohDiW7th_FuncMem ftohDiW7th_FuncTimer fun functional hasana hasanatomical hasatlas hasfreq hasfun hasmsk hasroi hastime i isUnstructuredFun max_op max_pow msk opacmax opacmin pow preamble_argin qi T tim tLabel tr TR trLabel 

clear global h tSlider trSlider ax counter

global h
global counter
counter = 1;
max_pow = 0;
min_pow = 0;
for i = 1:length(reldippow)

    if(max_pow < max(reldippow(i).pow))
        max_pow = max(reldippow(i).pow);
    end
    if(min_pow > min(reldippow(i).pow))
        min_pow = min(reldippow(i).pow);
    end
end

functional.pow     = reldippow(1).pow;
functional.inside  = inside;
functional.pos     = pos;
functional.tri     = tri;
cfg                = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
max_op             = max_pow/2;
min_op             = min_pow/2;
cfg.funcolorlim    = [min_pow max_pow];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [-10 10];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'Surface_NYC.mat'; % Cortical sheet from canonical MNI brain
cfg.surfdownsample = 10;                % downsample to speed up processing

tim                = mean(timebin,2);
% tr                 = length(alldippow(1).pow(:,1));
% trialnum           = 1:tr;








% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
% ft_nargin   = nargin;
% ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
% ft_preamble init
ft_preamble debug
ft_preamble loadvar functional anatomical
ft_preamble provenance functional anatomical

% the ft_abort variable is set to true or false in ft_preamble_init
% if ft_abort
%   return
% end

% this is not supported any more as of 26/10/2011
if ischar(functional)
  ft_error('please use cfg.inputfile instead of specifying the input variable as a sting');
end

% ensure that old and unsupported options are not being relied on by the end-user's script
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'location', 'interactive', 'auto'});
% instead of specifying cfg.coordsys, the user should specify the coordsys in the data
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'coordsys', 'inputcoord', 'inputcoordsys', 'coordinates'});
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2837
cfg = ft_checkconfig(cfg, 'renamed', {'viewdim', 'axisratio'});
cfg = ft_checkconfig(cfg, 'renamed', {'newfigure', 'figure'});

if isfield(cfg, 'atlas') && ~isempty(cfg.atlas)
  % the atlas lookup requires the specification of the coordsys
  functional     = ft_checkdata(functional, 'datatype', {'source', 'volume'}, 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');
else
  % check if the input functional is valid for this function, a coordsys is not directly needed
  functional     = ft_checkdata(functional, 'datatype', {'source', 'volume'}, 'feedback', 'yes', 'hasunit', 'yes');
end

% set the defaults for all methods
cfg.method        = ft_getopt(cfg, 'method',        'ortho');
cfg.funparameter  = ft_getopt(cfg, 'funparameter',  []);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter', []);
cfg.maskstyle     = ft_getopt(cfg, 'maskstyle',     'opacity');
cfg.downsample    = ft_getopt(cfg, 'downsample',    1);
cfg.title         = ft_getopt(cfg, 'title',         []);
cfg.figurename    = ft_getopt(cfg, 'figurename',    []);
cfg.atlas         = ft_getopt(cfg, 'atlas',         []);
cfg.marker        = ft_getopt(cfg, 'marker',        []);
cfg.markersize    = ft_getopt(cfg, 'markersize',    5);
cfg.markercolor   = ft_getopt(cfg, 'markercolor',   [1 1 1]);
cfg.colorbar      = ft_getopt(cfg, 'colorbar',      'yes');
cfg.colorbartext  = ft_getopt(cfg, 'colorbartext',  '');
cfg.voxelratio    = ft_getopt(cfg, 'voxelratio',    'data'); % display size of the voxel, 'data' or 'square'
cfg.axisratio     = ft_getopt(cfg, 'axisratio',     'data'); % size of the axes of the three orthoplots, 'square', 'voxel', or 'data'
cfg.visible       = ft_getopt(cfg, 'visible',       'on');
cfg.clim          = ft_getopt(cfg, 'clim',          [0 1]); % this is used to scale the orthoplot
cfg.intersectmesh = ft_getopt(cfg, 'intersectmesh');
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');

if ~isfield(cfg, 'anaparameter')
  if isfield(functional, 'anatomy')
    cfg.anaparameter = 'anatomy';
  else
    cfg.anaparameter = [];
  end
end

% set the common defaults for the functional data
cfg.funcolormap   = ft_getopt(cfg, 'funcolormap',   'auto');
cfg.funcolorlim   = ft_getopt(cfg, 'funcolorlim',   'auto');

% set the common defaults for the statistical data
cfg.opacitymap    = ft_getopt(cfg, 'opacitymap',    'auto');
cfg.opacitylim    = ft_getopt(cfg, 'opacitylim',    'auto');
cfg.roi           = ft_getopt(cfg, 'roi',           []);
cfg.maskstyle     = ft_getopt(cfg, 'maskstyle',     'opacity');

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, functional);
cfg.maskparameter = parameterselection(cfg.maskparameter, functional);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end

if isfield(functional, 'time') || isfield(functional, 'freq')
  % make a selection of the time and/or frequency dimension
  tmpcfg = keepfields(cfg, {'frequency', 'avgoverfreq', 'keepfreqdim', 'latency', 'avgovertime', 'keeptimedim', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  functional = ft_selectdata(tmpcfg, functional);
  % restore the provenance information
  [cfg, functional] = rollback_provenance(cfg, functional);
end

% the data can be passed as input argument or can be read from disk
hasanatomical = exist('anatomical', 'var');

if hasanatomical && ~strcmp(cfg.method, 'cloud') % cloud method should be able to take multiple surfaces and does not require interpolation
  % interpolate on the fly, this also does the downsampling if requested
  tmpcfg = keepfields(cfg, {'downsample', 'interpmethod', 'sphereradius', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  tmpcfg.parameter = cfg.funparameter;
  functional = ft_sourceinterpolate(tmpcfg, functional, anatomical);
  [cfg, functional] = rollback_provenance(cfg, functional);
  cfg.anaparameter = 'anatomy';
elseif ~hasanatomical && cfg.downsample~=1
  % optionally downsample the functional volume
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  tmpcfg.parameter = {cfg.funparameter, cfg.maskparameter, cfg.anaparameter};
  functional = ft_volumedownsample(tmpcfg, functional);
  [cfg, functional] = rollback_provenance(cfg, functional);
end

if isfield(functional, 'dim') && isfield(functional, 'transform')
  % this is a regular 3D functional volume
  isUnstructuredFun = false;

  % align the volume's coordinate system approximately to the voxels axes, this puts the box upright
  functional = align_ijk2xyz(functional);

elseif isfield(functional, 'dim') && isfield(functional, 'pos')
  % these are positions that can be mapped onto a 3D regular grid
  isUnstructuredFun  = false;
  % construct the transformation matrix from the positions
  functional.transform = pos2transform(functional.pos, functional.dim);
else
  % this is functional data on irregular positions, such as a cortical sheet
  isUnstructuredFun = true;
end

% this only relates to the dimensions of the geometry, which is npos*1 or nx*ny*nz
if isUnstructuredFun
  dim = [size(functional.pos,1) 1];
else
  dim = functional.dim;
end


% get the elements that will be plotted
hasatlas = ~isempty(cfg.atlas);
if hasatlas
  [atlas, functional] = handle_atlas_input(cfg.atlas, functional);
end

hasroi = ~isempty(cfg.roi);
if hasroi
  if ~hasatlas
    ft_error('specify cfg.atlas which specifies cfg.roi')
  else
    % get the mask
    tmpcfg = keepfields(cfg, {'roi', 'atlas'});
    roi = ft_volumelookup(tmpcfg, functional);
  end
end

hasana = isfield(functional, cfg.anaparameter);
if hasana
  ana = getsubfield(functional, cfg.anaparameter);
  if isa(ana, 'uint8') || isa(ana, 'uint16') || isa(ana, 'int8') || isa(ana, 'int16')
    ana = double(ana);
  end
  fprintf('scaling anatomy to [0 1]\n');
  dmin = min(ana(:));
  dmax = max(ana(:));
  ana  = (ana-dmin)./(dmax-dmin);
  ana  = reshape(ana, dim);
end

%%% funparameter
hasfun = isfield(functional, cfg.funparameter);
if hasfun
  fun = getsubfield(functional, cfg.funparameter);
    
  dimord = getdimord(functional, cfg.funparameter);
  dimtok = tokenize(dimord, '_');
  
  % replace the cell-array functional with a normal array
  if strcmp(dimtok{1}, '{pos}')
    tmpdim = getdimsiz(functional, cfg.funparameter);
    tmpfun = nan(tmpdim);
    insideindx = find(functional.inside);
    for i=insideindx(:)'
      tmpfun(i,:) = fun{i};
    end
    fun = tmpfun;
    clear tmpfun
    dimtok{1} = 'pos';  % update the description of the dimensions
    dimord([1 5]) = []; % remove the { and }
  end
  
  % ensure that the functional data is real
  if ~isreal(fun)
    ft_warning('functional data is complex, taking absolute value');
    fun = abs(fun);
  end
  
  if ~isa(fun, 'double')
    ft_warning('converting functional data to double precision');
    fun = double(fun);
  end
  
  if strcmp(dimord, 'pos_rgb') || (ndims(fun)>3 && size(fun,4)==3)
    % treat functional data as rgb values
    if any(fun(:)>1 | fun(:)<0)
      % scale
      tmpdim = size(fun);
      nvox   = prod(tmpdim(1:end-1));
      tmpfun = reshape(fun,[nvox tmpdim(end)]);
      m1     = max(tmpfun,[],1);
      m2     = min(tmpfun,[],1);
      tmpfun = (tmpfun-m2(ones(nvox,1),:))./(m1(ones(nvox,1),:)-m2(ones(nvox,1),:));
      fun    = reshape(tmpfun, tmpdim);
      clear tmpfun
    end
    qi      = 1;
    hasfreq = 0;
    hastime = 0;
    
    doimage = 1;
    fcolmin = 0;
    fcolmax = 1;
    
    cfg.funcolormap = 'rgb';
  else
    % determine scaling min and max (fcolmin fcolmax) and funcolormap
    if ~isa(fun, 'logical')
      funmin = min(fun(:));
      funmax = max(fun(:));
    else
      funmin = 0;
      funmax = 1;
    end
    % smart automatic limits
    if isequal(cfg.funcolorlim, 'auto')
      if sign(funmin)>-1 && sign(funmax)>-1
        cfg.funcolorlim = 'zeromax';
      elseif sign(funmin)<1 && sign(funmax)<1
        cfg.funcolorlim = 'minzero';
      else
        cfg.funcolorlim = 'maxabs';
      end
    end
    if ischar(cfg.funcolorlim)
      % limits are given as string
      if isequal(cfg.funcolorlim, 'maxabs')
        fcolmin = -max(abs([funmin,funmax]));
        fcolmax =  max(abs([funmin,funmax]));
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'default'; end
      elseif isequal(cfg.funcolorlim, 'zeromax')
        fcolmin = 0;
        fcolmax = funmax;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'hot'; end
      elseif isequal(cfg.funcolorlim, 'minzero')
        fcolmin = funmin;
        fcolmax = 0;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'cool'; end
      else
        ft_error('do not understand cfg.funcolorlim');
      end
    else
      % limits are numeric
      fcolmin = cfg.funcolorlim(1);
      fcolmax = cfg.funcolorlim(2);
      % smart colormap
      if isequal(cfg.funcolormap, 'auto')
        if sign(fcolmin) == -1 && sign(fcolmax) == 1
          cfg.funcolormap = 'default';
        else
          if fcolmin < 0
            cfg.funcolormap = 'cool';
          else
            cfg.funcolormap = 'hot';
          end
        end
      end
    end % if ischar
    clear funmin funmax
    
    % what if fun is 4D?
    if ndims(fun)>3 || prod(dim)==size(fun,1)
      if strcmp(dimord, 'pos_freq_time') || strcmp(dimord, 'dim1_dim2_dim3_freq_time')
        % functional contains time-frequency representation
        qi      = [1 1];
        hasfreq = numel(functional.freq)>1;
        hastime = numel(functional.time)>1;
        fun     = reshape(fun, [dim numel(functional.freq) numel(functional.time)]);
      elseif strcmp(dimord, 'pos_time') || strcmp(dimord, 'dim1_dim2_dim3_time')
        % functional contains evoked field
        qi      = 1;
        hasfreq = 0;
        hastime = numel(functional.time)>1;
        fun     = reshape(fun, [dim numel(functional.time)]);
      elseif strcmp(dimord, 'pos_ori_time') || strcmp(dimord, 'dim1_dim2_dim3_ori_time')
        % functional contains evoked field
        qi      = 1;
        hasfreq = 0;
        hastime = numel(functional.time)>1;
        % the following will fail if the number of orientations is larger than 1
        fun     = reshape(fun, [dim numel(functional.time)]);
      elseif strcmp(dimord, 'pos_freq') || strcmp(dimord, 'dim1_dim2_dim3_freq')
        % functional contains frequency spectra
        qi      = 1;
        hasfreq = numel(functional.freq)>1;
        hastime = 0;
        fun     = reshape(fun, [dim numel(functional.freq)]);
      else
        % functional contains scalar value for each position
        qi      = 1;
        hasfreq = 0;
        hastime = 0;
        fun     = reshape(fun, dim);
      end
    else
      % do nothing
      qi      = 1;
      hasfreq = 0;
      hastime = 0;
    end
    
    doimage = 0;
  end % if dimord has rgb or something else
  
else
  % there is no functional data
  qi      = 1;
  hasfreq = 0;
  hastime = 0;
  doimage = 0;
  fcolmin = 0; % needs to be defined for callback
  fcolmax = 1;
end

hasmsk = issubfield(functional, cfg.maskparameter);
if hasmsk
  if ~hasfun
    ft_error('you can not have a mask without functional data')
  else
    msk = getsubfield(functional, cfg.maskparameter);
    if islogical(msk) % otherwise sign() not posible
      msk = double(msk);
    end
  end
  % reshape to match fun
  if strcmp(dimord, 'pos_freq_time')
    % functional contains timefrequency representation
    msk     = reshape(msk, [dim numel(functional.freq) numel(functional.time)]);
  elseif strcmp(dimord, 'pos_time')
    % functional contains evoked field
    msk     = reshape(msk, [dim numel(functional.time)]);
  elseif strcmp(dimord, 'pos_freq')
    % functional contains frequency spectra
    msk     = reshape(msk, [dim numel(functional.freq)]);
  else
    msk     = reshape(msk, dim);
  end
  
  % determine scaling and opacitymap
  mskmin = min(msk(:));
  mskmax = max(msk(:));
  % determine the opacity limits and the opacity map
  % smart limits: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
  if isequal(cfg.opacitylim, 'auto')
    if isequal(cfg.funparameter,cfg.maskparameter)
      cfg.opacitylim = cfg.funcolorlim;
    else
      if sign(mskmin)>-1 && sign(mskmax)>-1
        cfg.opacitylim = 'zeromax';
      elseif sign(mskmin)<1 && sign(mskmax)<1
        cfg.opacitylim = 'minzero';
      else
        cfg.opacitylim = 'maxabs';
      end
    end
  end
  if ischar(cfg.opacitylim)
    % limits are given as string
    switch cfg.opacitylim
      case 'zeromax'
        opacmin = 0;
        opacmax = mskmax;
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'rampup'; end
      case 'minzero'
        opacmin = mskmin;
        opacmax = 0;
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'rampdown'; end
      case 'maxabs'
        opacmin = -max(abs([mskmin, mskmax]));
        opacmax =  max(abs([mskmin, mskmax]));
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'vdown'; end
      otherwise
        ft_error('incorrect specification of cfg.opacitylim');
    end % switch opacitylim
  else
    % limits are numeric
    opacmin = cfg.opacitylim(1);
    opacmax = cfg.opacitylim(2);
    if isequal(cfg.opacitymap, 'auto')
      if sign(opacmin)>-1 && sign(opacmax)>-1
        cfg.opacitymap = 'rampup';
      elseif sign(opacmin)<1 && sign(opacmax)<1
        cfg.opacitymap = 'rampdown';
      else
        cfg.opacitymap = 'vdown';
      end
    end
  end % handling opacitylim and opacitymap
  clear mskmin mskmax
else
  opacmin = [];
  opacmax = [];
end

% prevent outside fun from being plotted
if hasfun && ~hasmsk && isfield(functional, 'inside')
  hasmsk = 1;
  msk = zeros(dim);
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
  % make intelligent mask
  if isequal(cfg.method, 'surface')
    msk(functional.inside&isfinite(functional.(cfg.funparameter))) = 1;
    if any(functional.inside&~isfinite(functional.(cfg.funparameter)))
      ft_warning('functional data contains %d NaNs labeled as inside', sum(functional.inside&~isfinite(functional.(cfg.funparameter))));
    end
  else
    if hasana
      msk(functional.inside) = 0.5; % so anatomy is visible
    else
      msk(functional.inside) = 1;
    end
  end
end

% if region of interest is specified, mask everything besides roi
if hasfun && hasroi && ~hasmsk
  hasmsk = 1;
  msk = roi;
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
elseif hasfun && hasroi && hasmsk
  msk = roi .* msk;
  opacmin = [];
  opacmax = []; % has to be defined
elseif hasroi
  ft_error('you can not have a roi without functional data')
end

% %% give some feedback
% if ~hasfun && ~hasana
%   % this seems to be a problem that people often have due to incorrect specification of the cfg
%   ft_error('no anatomy is present and no functional data is selected, please check your cfg.funparameter');
% end
% if ~hasana
%   fprintf('not plotting anatomy\n');
% end
% if ~hasfun
%   fprintf('not plotting functional data\n');
% end
% if ~hasmsk
%   fprintf('not applying a mask on the functional data\n');
% end
% if ~hasatlas
%   fprintf('not using an atlas\n');
% end
% if ~hasroi
%   fprintf('not using a region-of-interest\n');
% end

% start building the figure

% open a new figure with the specified settings
h = open_figure_inourcode(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));
set(h, 'color', [1 1 1]);

global tSlider
tSlider = uislider(h,"Position",[100 125 600 3], "ValueChangingFcn",@(src,event)tSliderValueChanged(src,event,tim,reldippow));

tSlider.Limits = [min(tim), max(tim)];
tSlider.Value = min(tim);

tLabel = uilabel(h,"Position",[370 75 100 20], "Text", "Time (ms)" );

% global trSlider
% 
% trSlider = uislider(h,"Position",[100 130 600 3], "ValueChangingFcn",@(src,event)tSliderValueChanged(src,event,tim,trialnum,alldippow));
% trLabel = uilabel(h,"Position",[395 80 100 20], "Text", 'Trial (number)' );
% trSlider.Limits = [1 tr];
% trSlider.Value = 1;

% global meanbutton
% global bg
% bg = uibuttongroup(h,'Position',[295 20 120 40],'SelectionChangedFcn', @(src,event)tSliderValueChanged(src,event,tim,trialnum,alldippow));
% meanbutton = uitogglebutton(bg,'Position',[10 10 100 22]);
% meanbutton.Text = 'Average';
% meanbutton.Value = false;

global ax
ax = uiaxes(h,"Position",[150 175 500 500]);


% Plot the initial value

ft_sourceplot_update(tim, reldippow, cfg, dim, doimage, fcolmax, fcolmin, fun, functional, hasana, hasanatomical, hasatlas, hasfreq, hasfun, hasmsk, hastime, isUnstructuredFun, msk, opacmax, opacmin, qi)






%%

function tSliderValueChanged(src,event,tim,reldippow)
%     global h
%     clear global ax
    global ax
    global tSlider
%     global trSlider
%     global meanbutton
%     ax = uiaxes(h,"Position",[150 250 500 500]);

%     % hold off ax
%     ax.NextPlot = 'replace';


%     hold(ax,'off')
%     disp(tSlider.Value)
%         disp('here ')
%     ft_sourceplot_update(tim,trialnum,alldippow,cfg, dim, doimage, fcolmax, fcolmin, fun, functional, hasana, hasanatomical, hasatlas, hasfreq, hasfun, hasmsk, hastime, isUnstructuredFun, msk, opacmax, opacmin, qi)
    
    
    [~,I] = min(abs(tim-event.Value));
%     [~,IR] = min(abs(trialnum-trSlider.Value));
%     functional.pow     = alldippow(1,I).pow(IR,:);

% if (meanbutton == 1)
    ax.Children(2, 1).FaceVertexCData     = reldippow(I).pow';
    ax.Children(2, 1).FaceVertexAlphaData = reldippow(1,I).pow';
% elseif (meanbutton == 0)
%     ax.Children(2, 1).FaceVertexCData     = alldippow(1,I).pow(IR,:)';
%     ax.Children(2, 1).FaceVertexAlphaData = alldippow(1,I).pow(IR,:)';
% end
    
%     disp(['tSlider = ', num2str(event.Value)])
%         disp(['T = ', num2str(T)])
%         disp(['I = ', num2str(I)])
%         disp(['trSlider = ', num2str(trSlider.Value)])
% %         disp(['TR = ', num2str(TR)])
%         disp(['IR = ', num2str(IR)])

    drawnow;
end
%%

function ft_sourceplot_update(tim,reldippow, cfg, dim, doimage, fcolmax, fcolmin, fun, functional, hasana, hasanatomical, hasatlas, hasfreq, hasfun, hasmsk, hastime, isUnstructuredFun, msk, opacmax, opacmin, qi)
% ft_preamble init

% set color and opacity mapping for this figure
global h
global counter
global tSlider
% global trSlider
global ax
% global meanbutton

%     T=tSlider.Value;
%     TR=trSlider.Value;
   
%     [~,I] = min(abs(tim-tSlider.Value));
%     [~,IR] = min(abs(trialnum-trSlider.Value));
%     if (meanbutton == 1)
    functional.pow     = reldippow(1).pow;
%     elseif(meanbutton == 0)
%         functional.pow     = alldippow(1,I).pow(IR,:);
%     end
    
    %     if (counter == 2)
%         disp(['tSlider = ', num2str(tSlider.Value)])
%         disp(['T = ', num2str(T)])
%         disp(['I = ', num2str(I)])
%         disp(['trSlider = ', num2str(trSlider.Value)])
%         disp(['TR = ', num2str(TR)])
%         disp(['IR = ', num2str(IR)])
%         disp('here here')
%     end
    

if hasfun
  if ischar(cfg.funcolormap) && ~strcmp(cfg.funcolormap, 'rgb')
    cfg.funcolormap = ft_colormap_inourcode(cfg.funcolormap);
  elseif iscell(cfg.funcolormap)
    cfg.funcolormap = ft_colormap(cfg.funcolormap{:});
  end
end
if hasmsk
  cfg.opacitymap = alphamap_inourcode(cfg.opacitymap);
  alphamap_inourcode(cfg.opacitymap);
  if ndims(fun)>3 && ndims(msk)==3 && ~isequal(cfg.funcolormap, 'rgb')
    siz = size(fun);
    msk = repmat(msk, [1 1 1 siz(4:end)]);
  end
end

switch cfg.method
  case 'slice'
    assert(~hastime, 'method "%s" does not support time', cfg.method);
    assert(~hasfreq, 'method "%s" does not support freq', cfg.method);
    
    % set the defaults for method=slice
    cfg.nslices    = ft_getopt(cfg, 'nslices',    20);
    cfg.slicedim   = ft_getopt(cfg, 'slicedim',   3);
    cfg.slicerange = ft_getopt(cfg, 'slicerange', 'auto');
    
    % ADDED BY JM: allow for slicedim different than 3
    switch cfg.slicedim
      case 1
        if hasana, ana = permute(ana,[2 3 1]); end
        if hasfun, fun = permute(fun,[2 3 1]); end
        if hasmsk, msk = permute(msk,[2 3 1]); end
        cfg.slicedim=3;
      case 2
        if hasana, ana = permute(ana,[3 1 2]); end
        if hasfun, fun = permute(fun,[3 1 2]); end
        if hasmsk, msk = permute(msk,[3 1 2]); end
        cfg.slicedim=3;
      otherwise
        % nothing needed
    end
    
    %%%%% select slices
    if ~ischar(cfg.slicerange)
      ind_fslice = cfg.slicerange(1);
      ind_lslice = cfg.slicerange(2);
    elseif isequal(cfg.slicerange, 'auto')
      if hasfun % default
        if isfield(functional, 'inside')
          
          insideMask = false(size(fun));
          insideMask(functional.inside) = true;
          
          ind_fslice = find(max(max(insideMask,[],1),[],2), 1, 'first');
          ind_lslice = find(max(max(insideMask,[],1),[],2), 1, 'last');
        else
          ind_fslice = find(~isnan(max(max(fun,[],1),[],2)), 1, 'first');
          ind_lslice = find(~isnan(max(max(fun,[],1),[],2)), 1, 'last');
        end
      elseif hasana % if only ana, no fun
        ind_fslice = find(max(max(ana,[],1),[],2), 1, 'first');
        ind_lslice = find(max(max(ana,[],1),[],2), 1, 'last');
      else
        ft_error('no functional parameter and no anatomical parameter, can not plot');
      end
    else
      ft_error('do not understand cfg.slicerange');
    end
    ind_allslice = linspace(ind_fslice,ind_lslice,cfg.nslices);
    ind_allslice = round(ind_allslice);
    % make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
    if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end
    if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end
    if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end
    
    % update the dimensions of the volume
    if hasana
      dim=size(ana);
    else
      dim=size(fun);
    end
    
    %%%%% make a "quilt", that contain all slices on 2D patched sheet
    % Number of patches along sides of Quilt (M and N)
    % Size (in voxels) of side of patches of Quilt (m and n)
    
    % take care of a potential singleton 3rd dimension
    if numel(dim)<3
      dim(end+1:3) = 1;
    end
    
    m = dim(1);
    n = dim(2);
    M = ceil(sqrt(dim(3)));
    N = ceil(sqrt(dim(3)));
    
    num_patch = N*M;
    num_slice = (dim(cfg.slicedim));
    % put empty slides on ana, fun, msk, mskana to fill Quilt up
    if hasana; ana(:,:,end+1:num_patch)=0; end
    if hasfun; fun(:,:,end+1:num_patch)=0; end
    if hasmsk; msk(:,:,end+1:num_patch)=0; end
    % if hasmskana; mskana(:,:,end:num_patch)=0; end
    % put the slices in the quilt
    for iSlice = 1:num_slice
      xbeg = floor((iSlice-1)./M);
      ybeg = mod(iSlice-1, M);
      if hasana
        quilt_ana(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=ana(:,:,iSlice);
      end
      if hasfun
        quilt_fun(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=fun(:,:,iSlice);
      end
      if hasmsk
        quilt_msk(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=msk(:,:,iSlice);
      end
    end
    % make vols and scales, containes volumes to be plotted (fun, ana, msk), added by ingnie
    if hasana; vols2D{1} = quilt_ana; scales{1} = []; end % needed when only plotting ana
    if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end
    if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end
    
    % the transpose is needed for displaying the matrix using the MATLAB image() function
    if hasana;             ana = vols2D{1}'; end
    if hasfun && ~doimage; fun = vols2D{2}'; end
    if hasfun &&  doimage; fun = permute(vols2D{2},[2 1 3]); end
    if hasmsk;             msk = vols2D{3}'; end
    
    if hasana
      % scale anatomy between 0 and 1
      fprintf('scaling anatomy\n');
      amin = min(ana(:));
      amax = max(ana(:));
      ana = (ana-amin)./(amax-amin);
      clear amin amax;
      % convert anatomy into RGB values
      ana = cat(3, ana, ana, ana);
      ha = imagesc(ana);
    end
    hold on
    
    if hasfun
      
      if doimage
        hf = image(fun);
      else
        hf = imagesc(fun);
        try
          caxis(scales{2});
        end
        % apply the opacity mask to the functional data
        if hasmsk
          % set the opacity
          set(hf, 'AlphaData', msk)
          set(hf, 'AlphaDataMapping', 'scaled')
          try
            alim(scales{3});
          end
        elseif hasana
          set(hf, 'AlphaData', 0.5)
        end
        
      end
    end
    
    axis equal
    axis tight
    axis xy
    axis off
    
    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB colorbar
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
        ylabel(hc, cfg.colorbartext);
      else
        ft_warning('no colorbar possible without functional data')
      end
    end
    
  case 'ortho'
    % set the defaults for method=ortho
    cfg.location            = ft_getopt(cfg, 'location',            'auto');
    cfg.locationcoordinates = ft_getopt(cfg, 'locationcoordinates', 'head');
    cfg.crosshair           = ft_getopt(cfg, 'crosshair',           'yes');
    cfg.axis                = ft_getopt(cfg, 'axis',                'on');
    cfg.queryrange          = ft_getopt(cfg, 'queryrange',          3);
    
    if ~ischar(cfg.location)
      if strcmp(cfg.locationcoordinates, 'head')
        % convert the headcoordinates location into voxel coordinates
        loc = inv(functional.transform) * [cfg.location(:); 1];
        loc = round(loc(1:3));
      elseif strcmp(cfg.locationcoordinates, 'voxel')
        % the location is already in voxel coordinates
        loc = round(cfg.location(1:3));
      else
        ft_error('you should specify cfg.locationcoordinates');
      end
    else
      if isequal(cfg.location, 'auto')
        if hasfun
          if isequal(cfg.funcolorlim, 'maxabs')
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'zeromax')
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'minzero')
            loc = 'min';
          else % if numerical
            loc = 'max';
          end
        else
          loc = 'center';
        end
      else
        loc = cfg.location;
      end
    end
    
    % determine the initial intersection of the cursor (xi yi zi)
    if ischar(loc) && strcmp(loc, 'min')
      if isempty(cfg.funparameter)
        ft_error('cfg.location is min, but no functional parameter specified');
      end
      [dummy, minindx] = min(fun(:));
      [xi, yi, zi] = ind2sub(dim, minindx);
    elseif ischar(loc) && strcmp(loc, 'max')
      if isempty(cfg.funparameter)
        ft_error('cfg.location is max, but no functional parameter specified');
      end
      [dummy, maxindx] = max(fun(:));
      [xi, yi, zi] = ind2sub(dim, maxindx);
    elseif ischar(loc) && strcmp(loc, 'center')
      xi = round(dim(1)/2);
      yi = round(dim(2)/2);
      zi = round(dim(3)/2);
    elseif ~ischar(loc)
      % using nearest instead of round ensures that the position remains within the volume
      xi = nearest(1:dim(1), loc(1));
      yi = nearest(1:dim(2), loc(2));
      zi = nearest(1:dim(3), loc(3));
    end
    
    if numel(dim)<3
      ft_error('the input source structure cannot be reshaped into a volumetric 3D representation');
    end
    
    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));
    
    % axes settings
    if strcmp(cfg.axisratio, 'voxel')
      % determine the number of voxels to be plotted along each axis
      axlen1 = dim(1);
      axlen2 = dim(2);
      axlen3 = dim(3);
    elseif strcmp(cfg.axisratio, 'data')
      % determine the length of the edges along each axis
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      axlen1 = norm(cp_head(2,:)-cp_head(1,:));
      axlen2 = norm(cp_head(4,:)-cp_head(1,:));
      axlen3 = norm(cp_head(5,:)-cp_head(1,:));
    elseif strcmp(cfg.axisratio, 'square')
      % the length of the axes should be equal
      axlen1 = 1;
      axlen2 = 1;
      axlen3 = 1;
    end
    
    % this is the size reserved for subplot h1, h2 and h3
    h1size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h1size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h2size(1) = 0.82*axlen2/(axlen1 + axlen2);
    h2size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h3size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h3size(2) = 0.82*axlen2/(axlen2 + axlen3);
    
    if strcmp(cfg.voxelratio, 'square')
      voxlen1 = 1;
      voxlen2 = 1;
      voxlen3 = 1;
    elseif strcmp(cfg.voxelratio, 'data')
      % the size of the voxel is scaled with the data
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      voxlen1 = norm(cp_head(2,:)-cp_head(1,:))/norm(cp_voxel(2,:)-cp_voxel(1,:));
      voxlen2 = norm(cp_head(4,:)-cp_head(1,:))/norm(cp_voxel(4,:)-cp_voxel(1,:));
      voxlen3 = norm(cp_head(5,:)-cp_head(1,:))/norm(cp_voxel(5,:)-cp_voxel(1,:));
    end
    
    %% the figure is interactive, add callbacks
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_quit);
    
    %% create figure handles
    
    % axis handles will hold the anatomical functional if present, along with labels etc.
    h1 = axes('position',[0.06                0.06+0.06+h3size(2) h1size(1) h1size(2)]);
    h2 = axes('position',[0.06+0.06+h1size(1) 0.06+0.06+h3size(2) h2size(1) h2size(2)]);
    h3 = axes('position',[0.06                0.06                h3size(1) h3size(2)]);
    
    set(h1, 'Tag', 'ik', 'Visible', cfg.axis, 'XAxisLocation', 'top');
    set(h2, 'Tag', 'jk', 'Visible', cfg.axis, 'YAxisLocation', 'right'); % after rotating in ft_plot_ortho this becomes top
    set(h3, 'Tag', 'ij', 'Visible', cfg.axis);
    
    set(h1, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h2, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h3, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    
    % create structure to be passed to gui
    opt               = [];
    opt.dim           = dim;
    opt.ijk           = [xi yi zi];
    opt.h1size        = h1size;
    opt.h2size        = h2size;
    opt.h3size        = h3size;
    opt.handlesaxes   = [h1 h2 h3];
    opt.handlesfigure = h;
    opt.axis          = cfg.axis;
    if hasatlas
      opt.atlas = atlas;
    end
    if hasana && ~strcmp(cfg.maskstyle, 'colormix')
      opt.ana = ana;
    elseif hasana && strcmp(cfg.maskstyle, 'colormix')
      opt.background = ana;
    elseif ~hasana && ~strcmp(cfg.maskstyle, 'colormix')
      % nothing needed
    elseif ~hasana && strcmp(cfg.maskstyle, 'colormix')
      opt.background = zeros(size(fun));
    end
    if hasfun
      opt.fun = fun;
    end
    if hasmsk
      opt.msk = msk;
    end
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.usedim        = (isUnstructuredFun==false);
    opt.usepos        = (isUnstructuredFun==true);
    opt.hasatlas      = hasatlas;
    opt.hasfreq       = hasfreq;
    opt.hastime       = hastime;
    opt.hasmsk        = hasmsk;
    opt.hasfun        = hasfun;
    opt.hasana        = isfield(opt, 'ana');
    opt.hasbackground = isfield(opt, 'background');
    opt.qi            = qi;
    opt.tag           = 'ik';
    opt.functional    = functional;
    opt.fcolmin       = fcolmin;
    opt.fcolmax       = fcolmax;
    opt.opacmin       = opacmin;
    opt.opacmax       = opacmax;
    opt.clim          = cfg.clim; % contrast limits for the anatomy, see ft_volumenormalise
    opt.colorbar      = cfg.colorbar;
    opt.colorbartext  = cfg.colorbartext;
    opt.queryrange    = cfg.queryrange;
    opt.funcolormap   = cfg.funcolormap;
    opt.crosshair     = istrue(cfg.crosshair);
    if ~isempty(cfg.intersectmesh)
      % the data will be plotted in voxel space, so transform the meshes
      % accordingly, assuming the same coordinate system as the anatomical
      if ~isa(cfg.intersectmesh, 'cell')
        cfg.intersectmesh = {cfg.intersectmesh};
      end
      for m = 1:numel(cfg.intersectmesh)
        opt.intersectmesh{m} = ft_transform_geometry(inv(functional.transform), cfg.intersectmesh{m});
      end
    end
    
    %% do the actual plotting
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
    fprintf('\n');
    fprintf('click left mouse button to reposition the cursor\n');
    fprintf('click and hold right mouse button to update the position while moving the mouse\n');
    fprintf('use the arrowkeys to navigate in the current axis\n');
    
    
  case 'surface'
    assert(~hastime, 'method "%s" does not support time', cfg.method);
    assert(~hasfreq, 'method "%s" does not support freq', cfg.method);
    
    % set the defaults for method=surface
    cfg.downsample     = ft_getopt(cfg, 'downsample',     1);
    cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
    cfg.surffile       = ft_getopt(cfg, 'surffile', 'surface_white_both.mat'); % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
    cfg.surfinflated   = ft_getopt(cfg, 'surfinflated',  []);
    cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
    cfg.projmethod     = ft_getopt(cfg, 'projmethod',    'nearest');
    cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
    cfg.camlight       = ft_getopt(cfg, 'camlight',      'yes');
    cfg.facecolor      = ft_getopt(cfg, 'facecolor',    []);
    cfg.vertexcolor    = ft_getopt(cfg, 'vertexcolor',   'curv'); % curvature-dependent mix of cortex_light and cortex_dark
    cfg.edgecolor      = ft_getopt(cfg, 'edgecolor',     'none');
    
    % determine whether the source functional already contains a triangulation
    interpolate2surf = 0;
    if ~isUnstructuredFun
      % no triangulation present: interpolation should be performed
      fprintf('The source functional is defined on a 3D grid, interpolation to a surface mesh will be performed\n');
      interpolate2surf = 1;
    elseif isUnstructuredFun && isfield(functional, 'tri')
      fprintf('The source functional is defined on a triangulated surface, using the surface mesh description in the functional\n');
    elseif isUnstructuredFun
      % add a transform field to the functional
      fprintf('The source functional does not contain a triangulated surface, we may need to interpolate to a surface mesh\n');
      functional.transform = pos2transform(functional.pos);
      interpolate2surf = 1;
    end
    
    if interpolate2surf
      % deal with the interpolation
      % FIXME this should be dealt with by ft_sourceinterpolate
      
      % read the triangulated cortical surface from file
      surf = ft_read_headshape(cfg.surffile);
      if isfield(surf, 'transform')
        % compute the vertex positions in head coordinates
        surf.pos = ft_warp_apply(surf.transform, surf.pos);
      end
      % ensure that the surface is formatted properly
      surf = ft_checkdata(surf, 'hasunit', isfield(functional, 'unit'), 'hascoordsys', isfield(functional, 'coordsys'));
      if isfield(functional, 'unit')
        % ensure that the units are consistent, convert the units if required
        surf = ft_convert_units(surf, functional.unit);
      end
      if isfield(functional, 'coordsys') && isfield(surf, 'coordsys')
        % ensure that the coordinate systems match
        functional = fixcoordsys(functional);
        surf       = fixcoordsys(surf);
        assert(isequal(functional.coordsys, surf.coordsys), 'coordinate systems do not match');
      else
        ft_notice('assuming that the coordinate systems match');
      end
      
      % downsample the cortical surface
      if cfg.surfdownsample > 1
        if ~isempty(cfg.surfinflated)
          ft_error('downsampling the surface is not possible in combination with an inflated surface');
        end
        fprintf('downsampling surface from %d vertices\n', size(surf.pos,1));
        [temp.tri, temp.pos] = reducepatch(surf.tri, surf.pos, 1/cfg.surfdownsample);
        % find indices of retained patch faces
        [dummy, idx] = ismember(temp.pos, surf.pos, 'rows');
        idx(idx==0)  = [];
        surf.tri = temp.tri;
        surf.pos = temp.pos;
        clear temp
        % downsample other fields
        if isfield(surf, 'curv'),       surf.curv       = surf.curv(idx);       end
        if isfield(surf, 'sulc'),       surf.sulc       = surf.sulc(idx);       end
        if isfield(surf, 'hemisphere'), surf.hemisphere = surf.hemisphere(idx); end
        if isfield(surf, 'inside'),     surf.inside     = surf.inside(idx);     end
      end
      
      % these are required
      if ~isfield(functional, 'inside')
        functional.inside = true(dim);
      end
      
      fprintf('%d voxels in functional data\n', prod(dim));
      fprintf('%d vertices in cortical surface\n', size(surf.pos,1));
      
      tmpcfg = [];
      tmpcfg.parameter = {cfg.funparameter};
      if ~isempty(cfg.maskparameter)
        % it was specified by the user
        tmpcfg.parameter = [tmpcfg.parameter {cfg.maskparameter}];
        maskparameter    = cfg.maskparameter;
      elseif hasmsk
        % it was constructed on the fly
        functional.mask  = msk;
        tmpcfg.parameter = [tmpcfg.parameter {'mask'}];
        maskparameter    = 'mask'; % temporary variable
      end
      tmpcfg.interpmethod = cfg.projmethod;
      tmpcfg.distmat      = cfg.distmat;
      tmpcfg.sphereradius = cfg.sphereradius;
      tmpcfg.projvec      = cfg.projvec;
      tmpcfg.projcomb     = cfg.projcomb;
      tmpcfg.projweight   = cfg.projweight;
      tmpcfg.projthresh   = cfg.projthresh;
      tmpdata             = ft_sourceinterpolate(tmpcfg, functional, surf);
      
      if hasfun, val      = getsubfield(tmpdata, cfg.funparameter);  val     = val(:);     end
      if hasmsk, maskval  = getsubfield(tmpdata, maskparameter);     maskval = maskval(:); end
      
      if ~isempty(cfg.projthresh) && hasmsk
        maskval(abs(val) < cfg.projthresh*max(abs(val(:)))) = 0;
      end
      
    else
      surf     = [];
      surf.pos = functional.pos;
      surf.tri = functional.tri;
      
      % if hasfun, val     = fun(functional.inside(:)); end
      % if hasmsk, maskval = msk(functional.inside(:)); end
      if hasfun, val     = fun(:); end
      if hasmsk, maskval = msk(:); end
      
    end
    
    if ~isempty(cfg.surfinflated)
      if ~isstruct(cfg.surfinflated)
        % read the inflated triangulated cortical surface from file
        surf = ft_read_headshape(cfg.surfinflated);
      else
        surf = cfg.surfinflated;
        if isfield(surf, 'transform')
          % compute the surface vertices in head coordinates
          surf.pos = ft_warp_apply(surf.transform, surf.pos);
        end
      end
    end
    
    %------do the plotting
    if ~hasfun
      ft_plot_mesh_inourcode(surf,'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
      
    elseif hasfun
      if ~hasmsk || all(maskval(:)==1)
        ft_plot_mesh_inourcode(surf,'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val, 'clim', [fcolmin fcolmax], 'colormap', cfg.funcolormap);
        
      elseif hasmsk
        switch cfg.maskstyle
          case 'opacity'
            % this needs to be handled here, rather than in ft_plot_mesh,
            % because the latter function needs to be called twice: once
            % for drawing the background, with possibly a user-specified
            % background color, and the second time for the functional data
            ft_plot_mesh_inourcode(surf, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
            ft_plot_mesh_inourcode(surf, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val, 'facealpha', maskval, 'clim', [fcolmin fcolmax], 'alphalim', [opacmin opacmax], 'alphamap_inourcode', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
            
          case 'colormix'
            % convert the specification of the background color + functional
            % color + opacity into a single rgb value to speed up the rendering
            ft_plot_mesh_inourcode(surf, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val, 'facealpha', maskval, 'clim', [fcolmin fcolmax], 'alphalim', [opacmin opacmax], 'alphamap_inourcode', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'colormix');
            
        end
      end
    end
    
%     if (counter == 1)
    
        lighting(ax,'gouraud')
%     end
if (counter == 1)
          
          counter = 2;
      end
    if istrue(cfg.camlight)
      camlight(ax)
      
    end
    
    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB colorbar
        hc = colorbar(ax);
        if strcmp(cfg.maskstyle, 'opacity')
          % functional values are according to original input values
          set(hc, 'YLim', [fcolmin fcolmax]);
          ylabel(hc, cfg.colorbartext);
        else
          % functional values have been transformed to be scaled
        end
      else
        ft_warning('no colorbar possible without functional data')
      end
    end
    
  case 'glassbrain'
    assert(~hastime, 'method "%s" does not support time', cfg.method);
    assert(~hasfreq, 'method "%s" does not support freq', cfg.method);
    
    % This is implemented using a recursive call with an updated functional data
    % structure. The functional volume is replaced by a volume in which the maxima
    % are projected to the "edge" of the volume.
    
    tmpfunctional = keepfields(functional, {'dim', 'transform'});
    
    if hasfun
      if isfield(functional, 'inside')
        fun(~functional.inside) = nan;
      end
      fun(1,:,:) = max(fun, [], 1); % get the projection along the 1st dimension
      fun(:,1,:) = max(fun, [], 2); % get the projection along the 2nd dimension
      fun(:,:,1) = max(fun, [], 3); % get the projection along the 3rd dimension
      tmpfunctional.(cfg.funparameter) = fun;
    end
    
    if hasana
      if isfield(functional, 'inside')
        % ana(~functional.inside) = nan;
      end
      ana(1,:,:) = max(ana, [], 1); % get the projection along the 1st dimension
      ana(:,1,:) = max(ana, [], 2); % get the projection along the 2nd dimension
      ana(:,:,1) = max(ana, [], 3); % get the projection along the 3rd dimension
      tmpfunctional.(cfg.anaparameter) = ana;
    end
    
    tmpcfg                      = keepfields(cfg, {'anaparameter', 'funparameter', 'funcolorlim', 'funcolormap', 'opacitylim', 'axis', 'renderer', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
    tmpcfg.method               = 'ortho';
    tmpcfg.location             = [1 1 1];
    tmpcfg.locationcoordinates  = 'voxel';
    ft_sourceplot(tmpcfg, tmpfunctional);
    
  case 'vertex'
    assert(~hastime, 'method "%s" does not support time', cfg.method);
    assert(~hasfreq, 'method "%s" does not support freq', cfg.method);
    
    if isUnstructuredFun
      pos = functional.pos;
    else
      [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      pos = ft_warp_apply(functional.transform, [X(:) Y(:) Z(:)]);
    end
    
    if isfield(functional, 'inside')
      pos = pos(functional.inside,:);
      if hasfun
        fun = fun(functional.inside);
      end
    end
    
    % scale the functional data between -30 and 30
    fun = 30*fun/max(abs(fun(:)));
    if any(fun<=0)
      ft_warning('using red for positive and blue for negative functional values')
      col = zeros(numel(fun), 3); % RGB
      col(fun>0,1) = 1;  % red
      col(fun<0,3) = 1;  % blue
      fun(fun==0) = eps; % these will be black
      ft_plot_mesh_inourcode(pos, 'vertexsize', abs(fun), 'vertexcolor', col);
    else
      ft_plot_mesh_inourcode(pos, 'vertexsize', fun, 'vertexcolor', 'k');
    end
    
    % ensure that the axes don't change if you rotate
    axis vis3d
    
  case 'cloud'
    assert(~hastime, 'method "%s" does not support time', cfg.method);
    assert(~hasfreq, 'method "%s" does not support freq', cfg.method);
    
    % some defaults depend on the geometrical units
    scale = ft_scalingfactor('mm', functional.unit);
    
    % set the defaults for method=cloud
    cfg.cloudtype          = ft_getopt(cfg, 'cloudtype', 'cloud');
    cfg.scalerad           = ft_getopt(cfg, 'scalerad', 'yes');
    cfg.ptsize             = ft_getopt(cfg, 'ptsize', 1);
    cfg.ptdensity          = ft_getopt(cfg, 'ptdensity', 20);
    cfg.ptgradient         = ft_getopt(cfg, 'ptgradient', .5);
    cfg.colorgrad          = ft_getopt(cfg, 'colorgrad', 'white');
    cfg.marker             = ft_getopt(cfg, 'marker', '.');
    cfg.slice              = ft_getopt(cfg, 'slice', 'none');
    cfg.ori                = ft_getopt(cfg, 'ori', 'y');
    cfg.slicepos           = ft_getopt(cfg, 'slicepos', 'auto');
    cfg.nslices            = ft_getopt(cfg, 'nslices', 1);
    cfg.minspace           = ft_getopt(cfg, 'minspace', 1);
    cfg.intersectcolor     = ft_getopt(cfg, 'intersectcolor', {'k'});
    cfg.intersectlinestyle = ft_getopt(cfg, 'intersectlinestyle', {'-'});
    cfg.intersectlinewidth = ft_getopt(cfg, 'intersectlinewidth', 2);
    cfg.ncirc              = ft_getopt(cfg, 'ncirc', 15);
    cfg.scalealpha         = ft_getopt(cfg, 'scalealpha', 'no');
    cfg.facecolor          = ft_getopt(cfg, 'facecolor', [0.781 0.762 0.664]);
    cfg.edgecolor          = ft_getopt(cfg, 'edgecolor', 'none');
    cfg.facealpha          = ft_getopt(cfg, 'facealpha', 1);
    cfg.edgealpha          = ft_getopt(cfg, 'edgealpha', 0);
    cfg.vertexcolor        = ft_getopt(cfg, 'vertexcolor', 'curv'); % curvature-dependent mix of cortex_light and cortex_dark
    if ~hasanatomical; anatomical = {}; end
    
    if isUnstructuredFun
      pos = functional.pos;
    else
      [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      pos = ft_warp_apply(functional.transform, [X(:) Y(:) Z(:)]);
    end
    
    if hasmsk
      pos = pos(logical(msk),:);
      if hasfun
        fun = fun(logical(msk));
      end
    end
    
    if strcmp(cfg.cloudtype, 'cloud') || strcmp(cfg.cloudtype, 'surf')
      % set the defaults for cloudtype=cloud & cloudtype=surf
      cfg.radius = ft_getopt(cfg, 'radius', 4*scale);
      cfg.rmin   = ft_getopt(cfg, 'rmin', 1*scale);
      
      ft_plot_cloud(pos, fun, 'mesh', anatomical,...
        'radius', cfg.radius, 'rmin', cfg.rmin, 'scalerad', cfg.scalerad, ...
        'ptsize', cfg.ptsize, 'ptdensity', cfg.ptdensity, 'ptgradient', cfg.ptgradient,...
        'colorgrad', cfg.colorgrad, 'colormap', cfg.funcolormap, 'clim', [fcolmin fcolmax], ...
        'unit', functional.unit, 'slice', cfg.slice, 'cloudtype', cfg.cloudtype, ...
        'ori', cfg.ori, 'slicepos', cfg.slicepos, 'nslices', cfg.nslices, 'minspace', cfg.minspace,...
        'intersectcolor', cfg.intersectcolor, 'intersectlinestyle', cfg.intersectlinestyle, ...
        'intersectlinewidth', cfg.intersectlinewidth, 'ncirc', cfg.ncirc, ...
        'scalealpha', cfg.scalealpha, 'facecolor', cfg.facecolor, 'edgecolor', cfg.edgecolor,...
        'facealpha', cfg.facealpha, 'edgealpha', cfg.edgealpha, 'marker', cfg.marker,...
        'vertexcolor', cfg.vertexcolor);
      
    elseif strcmp(cfg.cloudtype, 'point')
      if strcmp(cfg.slice, '2d') || strcmp(cfg.slice, '3d')
        error('slices are not supported for cloudtype=''point''')
      end
      
      % set the defaults for cloudtype=point
      cfg.radius = ft_getopt(cfg, 'radius', 40*scale);
      cfg.rmin   = ft_getopt(cfg, 'rmin', 10*scale);
      
      % functional scaling
      cmap    = cfg.funcolormap;
      cmid    = size(cmap,1)/2;                            % colorbar middle
      clim    = [fcolmin fcolmax];                         % color limits
      colscf  = 2*( (fun-clim(1)) / (clim(2)-clim(1)) )-1; % color between -1 and 1, bottom vs. top colorbar
      colscf(colscf>1)=1; colscf(colscf<-1)=-1;            % clip values outside the [-1 1] range
      radscf  = fun-(min(abs(fun)) * sign(max(fun)));      % radius between 0 and 1, small vs. large pos/neg effect
      radscf  = abs( radscf / max(abs(radscf)) );
      
      if strcmp(cfg.scalerad, 'yes')
        rmax = cfg.rmin+(cfg.radius-cfg.rmin)*radscf; % maximum radius of the clouds
      else
        rmax = ones(length(pos), 1)*cfg.radius; % each cloud has the same radius
      end
      
      % plot functional
      for n = 1:size(pos,1) % sensor loop
        indx  = ceil(cmid) + sign(colscf(n))*floor(abs(colscf(n)*cmid));
        indx  = max(min(indx,size(cmap,1)),1);  % index should fall within the colormap
        fcol  = cmap(indx,:);                   % color [Nx3]
        hold on; plot3(pos(n,1), pos(n,2), pos(n,3), 'Marker', cfg.marker, 'MarkerSize', rmax(n), 'Color', fcol, 'Linestyle', 'none');
      end
      
      % plot anatomical
      if hasanatomical
        ft_plot_mesh_inourcode(anatomical, 'facecolor', cfg.facecolor, 'EdgeColor', cfg.edgecolor, 'facealpha', cfg.facealpha, 'edgealpha', cfg.edgealpha, 'vertexcolor', cfg.vertexcolor);
        material dull
      end
      
      % color settings
      ft_colormap(cmap);
      if ~isempty(clim) && clim(2)>clim(1)
        caxis(h, clim);
      end
    end
    
    if istrue(cfg.colorbar)
      if ~strcmp(cfg.slice, '2d')
        c = colorbar;
      else % position the colorbar so that it does not change the axis of the last subplot
        subplotpos = get(subplot(cfg.nslices,1,cfg.nslices), 'Position'); % position of the bottom or rightmost subplot
        c = colorbar('Position', [subplotpos(1)+subplotpos(3)+0.01 subplotpos(2) .03 subplotpos(2)+subplotpos(4)*(cfg.nslices+.1)]);
      end
      ylabel(c, cfg.colorbartext);
    end
    
    
  otherwise
    ft_error('unsupported method "%s"', cfg.method);
end

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

% set the figure window title
if ~isempty(dataname)
  set(h, 'Name', sprintf('%d: %s: %s', double(h), mfilename, join_str(', ', dataname)));
else
  set(h, 'Name', sprintf('%d: %s', double(h), mfilename));
end
set(h, 'NumberTitle', 'off');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous functional
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, the subplots are well-controlled in this case
menu_fieldtrip(h, cfg, true);

% if ~ft_nargout
%   % don't return anything
%   clear cfg
% end

end
%%

function cmap = ft_colormap_inourcode(varargin)
global h
% FT_COLORMAP is a wrapper function with the same usage as the normal COLORMAP
% function, but it also knows about the colormaps from BREWERMAP and some colormaps
% from MATPLOTLIB. The recommended colormaps include 'parula', 'cividis', 'balance',
% and '*RdBu'.
%
% Use as
%   ft_colormap(name)
%   ft_colormap(name, n)
%   ft_colormap(handle, name)
%   ft_colormap(handle, name, n)
%
% The name is a string that specifies the colormap (see below). The optional handle
% can be used to specify the current figure (which is the default, see GCF) or the
% current axes (see GCA). The optional parameter n determines the number of steps or
% unique colors in the map (by default 64).
%
% The colormaps from MATLAB include 'parula', 'jet', 'hsv', 'hot', 'cool', 'spring',
% 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines',
% 'colorcube', 'prism', and 'flag'.
%
% The colormaps from MATPLOTLIB include 'cividis', 'inferno', 'magma', 'plasma',
% 'tab10', 'tab20', 'tab20b', 'tab20c', 'twilight', and 'viridis'.
%
% The colormaps from BREWERMAP include 'BrBG', 'PRGn', 'PiYG', 'PuOr', 'RdBu',
% 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1',
% 'Pastel2', 'Set1', 'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens',
% 'Greys', 'OrRd', 'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
% 'YlGn', 'YlGnBu', 'YlOrBr', and 'YlOrRd', plus their reverse when prefixed with '*'.
%
% The colormaps from CMOCEAN include 'thermal', 'haline', 'solar', 'ice', 'gray',
% 'oxy', 'deep', 'dense', 'algae', 'matter', 'turbid', 'speed', 'amp', 'tempo',
% 'rain', 'phase', 'topo', 'balance', 'delta', 'curl', 'diff', and 'tarn'.
%
% The colormaps from COLORCET include 'blueternary', 'coolwarm', 'cyclicgrey',
% 'depth', 'divbjy', 'fire', 'geographic', 'geographic2', 'gouldian', 'gray',
% 'greenternary', 'grey', 'heat', 'phase2', 'phase4', 'rainbow', 'rainbow2',
% 'rainbow3', 'rainbow4', 'redternary', 'reducedgrey', 'yellowheat', and all the ones
% with symbolic names.
%
% To reverse any of these these colormaps you can add a minus sign in front, like
% '-phase', '-balance' or '-RdBu'.
%
% Relevant publications:
% - Crameri et al. 2020. The misuse of colour in science communication. https://doi.org/10.1038/s41467-020-19160-7
% - Cooper et al. 2021. Over the rainbow: Guidelines for meaningful use of colour maps in neurophysiology. https://doi.org/10.1016/j.neuroimage.2021.118628
% - Kovesi 2015, Good colour maps: How to design them. https://doi.org/10.48550/arXiv.1509.03700
%
% See also COLORMAP, COLORMAPEDITOR, BREWERMAP, MATPLOTLIB, CMOCEAN, COLORCET

% Copyright (C) 2022, Jan-Mathijs Schoffelen, Robert Oostenveld
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

ft_hastoolbox('brewermap', 2);  % check and give a warning if it cannot be added
ft_hastoolbox('matplotlib', 2); % check and give a warning if it cannot be added
ft_hastoolbox('cmocean', 2);
ft_hastoolbox('colorcet', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpret the input arguments

% please note that ishandle('parula') returns a logical array
% hence the all(ishandle(...)), which in that case returns false

if nargin==0
  % apply the default colormap
  name   = 'default';
  n      = 64;

elseif nargin==1 && isnumeric(varargin{1})
  % the user specified an Nx3 array as colormap

  name   = varargin{1}; % note that this is not a string, it is dealt with below
  n      = nan;

elseif nargin==1 && all(ishandle(varargin{1}))
  % apply the default colormap on the specified handle
  handle = h;
  name   = 'default';
  n      = 64;

elseif nargin==1 && ischar(varargin{1})
  % apply the specified colormap on the current figure
  handle = h;
  name   = 'jet';
  n      = 64;

elseif nargin==2 && all(ishandle(varargin{1}))
  % apply the specified colormap on the specified handle
  handle = h;
  name   = varargin{2};
  n      = 64;

elseif nargin==2 && ischar(varargin{1})
  % apply the specified colormap with specified N on the current figure
  handle = h;
  name   = varargin{1};
  n      = varargin{2};

elseif nargin==3 && all(ishandle(varargin{1}))
  % apply the specified colormap with specified N on the specified handle
  handle = h;
  name   = varargin{2};
  n      = varargin{3};

else
  ft_error('incorrect input arguments');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the colormap

% the ones from brewermap are not m-files on disk
brewerlist = repmat(brewermap('list'), [2 1]);
% also include the reverse ones
for k = 1:numel(brewerlist)/2
  brewerlist{k} = sprintf('*%s', brewerlist{k});
end

% the ones from cmocean are not m-files on disk
cmoceanlist = {
  'thermal'
  'haline'
  'solar'
  'ice'
  'gray'
  'oxy'
  'deep'
  'dense'
  'algae'
  'matter'
  'turbid'
  'speed'
  'amp'
  'tempo'
  'rain'
  'phase'
  'topo'
  'balance'
  'delta'
  'curl'
  'diff'
  'tarn'
  };
% also include the reverse ones
for k = 1:numel(cmoceanlist)
  cmoceanlist{end+1} = sprintf('-%s', cmoceanlist{k});
end

% the ones from colorcet are not m-files on disk
colorcetlist = {
  'blueternary'
  'coolwarm'
  'cyclicgrey'
  'depth'
  'divbjy'
  'fire'
  'geographic'
  'geographic2'
  'gouldian'
  'gray'
  'greenternary'
  'grey'
  'heat'
  'phase2'
  'phase4'
  'rainbow'
  'rainbow2'
  'rainbow3'
  'rainbow4'
  'redternary'
  'reducedgrey'
  'yellowheat'
  };
% also include the reverse ones
for k = 1:numel(colorcetlist)
  colorcetlist{end+1} = sprintf('-%s', colorcetlist{k});
end

if isnumeric(name)
  % the user specified an Nx3 array as colormap
  cmap = name;
elseif ismember(name, brewerlist)
  cmap = brewermap(n, name);
elseif startsWith(name, '-') && ismember(name(2:end), brewerlist)
  cmap = brewermap(n, strrep(name, '-', '*'));
elseif ismember(name, cmoceanlist)
  cmap = cmocean(name, n);
elseif ismember(name, colorcetlist)
  cmap = colorcet(name, 'n', n, 'reverse', false);
elseif startsWith(name, '-') && ismember(name(2:end), colorcetlist)
  cmap = colorcet(name(2:end), 'n', n, 'reverse', true);
elseif isequal(name, 'default')
  % requires separate handling, because there's no function called default,
  % the default is taken care of by colormap
  cmap = name;
else
  % this works both for the MATLAB and the MATPLOTLIB colormaps
  % which have the different colormaps available as an m-file
  if name(1)=='-'
    cmap = flipud(feval(name(2:end), n));
  else
    cmap = feval(name, n);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the colormap and/or return it

if nargout>0
  % apply it to the specified figure or axis and return it as a Nx3 array
  colormap(h, cmap)
else
  % only apply it to the current figure
  colormap(h, cmap)
  clear cmap
end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
tag = get(curr_ax, 'tag');

functional = opt.functional;

h1 = opt.handlesaxes(1);
h2 = opt.handlesaxes(2);
h3 = opt.handlesaxes(3);

xi = opt.ijk(1);
yi = opt.ijk(2);
zi = opt.ijk(3);
qi = opt.qi;

if any([xi yi zi] > functional.dim) || any([xi yi zi] <= 0)
  return;
end

opt.ijk = [xi yi zi 1]';
if opt.usedim
  xyz = functional.transform * opt.ijk;
elseif opt.usepos
  ix  = sub2ind(opt.dim,xi,yi,zi);
  xyz = functional.pos(ix,:);
end
opt.ijk = opt.ijk(1:3);

% construct a string with user feedback
str1 = sprintf('voxel %d\nindices [%d %d %d]', sub2ind(functional.dim(1:3), xi, yi, zi), opt.ijk);

if isfield(functional, 'coordsys')
  cstr = sprintf('%s coordinates', functional.coordsys);
  [dirijk(1,:), dirijk(2,:), dirijk(3,:)] = coordsys2label(functional.coordsys, 1, 1);
  showcoordsys = ~isempty(functional.coordsys) && ~isequal(functional.coordsys, 'unknown');
else
  cstr = 'location';
  showcoordsys = false;
end
if isfield(functional, 'unit')
  switch functional.unit
    case 'm'
      ustr = sprintf('[%.3f %.3f %.3f] m', xyz(1:3));
    case 'cm'
      ustr = sprintf('[%.1f %.1f %.1f] cm', xyz(1:3));
    case 'mm'
      ustr = sprintf('[%.0f %.0f %.0f] mm', xyz(1:3));
    otherwise
      ustr = sprintf('[%f %f %f] %s', xyz(1:3), functional.unit);
  end
else
  ustr = sprintf('[%.3f %.3f %.3f]', xyz(1:3));
end
str2 = sprintf('%s %s', cstr, ustr);

if opt.hasfreq && opt.hastime
  str3 = sprintf('%.1f s, %.1f Hz', functional.time(opt.qi(2)), functional.freq(opt.qi(1)));
elseif ~opt.hasfreq && opt.hastime
  str3 = sprintf('%.1f s', functional.time(opt.qi(1)));
elseif opt.hasfreq && ~opt.hastime
  str3 = sprintf('%.1f Hz', functional.freq(opt.qi(1)));
else
  str3 = '';
end

if opt.hasfun
  if ~opt.hasfreq && ~opt.hastime
    val = opt.fun(xi, yi, zi);
  elseif ~opt.hasfreq && opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi);
  elseif opt.hasfreq && ~opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi);
  elseif opt.hasfreq && opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi(1), opt.qi(2));
  end
  str4 = sprintf('value %f', val);
else
  str4 = '';
end

if opt.hasatlas
  % determine the anatomical label of the current position
  lab = atlas_lookup(opt.atlas, (xyz(1:3)), 'coordsys', functional.coordsys, 'queryrange', opt.queryrange);
  if isempty(lab)
    lab = 'NA';
  else
    lab = unique(lab);
    tmp = sprintf('%s', strrep(lab{1}, '_', ' '));
    for i=2:length(lab)
      tmp = [tmp sprintf(', %s', strrep(lab{i}, '_', ' '))];
    end
    lab = tmp;
  end
else
  lab = 'NA';
end

if opt.hasana
  options = {'transform', eye(4),     'location', opt.ijk, 'style', 'subplot',...
             'update',    opt.update, 'doscale',  false,   'clim',  opt.clim};
  if isfield(opt, 'intersectmesh')
    options = cat(2, options, 'intersectmesh', opt.intersectmesh);
  end
  
  if opt.init
    tmph  = [h1 h2 h3];
    options = cat(2, options, {'parents', tmph});
    ft_plot_ortho(opt.ana, options{:});
    
    opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
    for i=1:length(opt.anahandles)
      opt.parenttag{i} = get(get(opt.anahandles(i), 'parent'), 'tag');
    end
    [i1,i2,i3] = intersect(opt.parenttag, {'ik' 'jk' 'ij'});
    opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
    opt.anahandles = opt.anahandles(:)';
    set(opt.anahandles, 'tag', 'ana');
    if isfield(opt, 'intersectmesh')
      opt.patchhandles = findobj(opt.handlesfigure, 'type', 'patch');
      opt.patchhandles = opt.patchhandles(i3(i2));
      opt.patchhandles = opt.patchhandles(:)';
      set(opt.patchhandles, 'tag', 'patch');
    end
  else
    options = cat(2, options, {'surfhandle', opt.anahandles});
    if isfield(opt, 'intersectmesh')
      options = cat(2, options, {'patchhandle', opt.patchhandles});
    end
    ft_plot_ortho(opt.ana, options{:});
  end
end

if opt.hasfun
  if opt.init
    tmph  = [h1 h2 h3];
    if isequal(opt.funcolormap, 'rgb')
      tmpfun = opt.fun;
      if opt.hasmsk
        tmpmask = opt.msk;
      end
    else
      tmpqi  = [opt.qi 1];
      tmpfun = opt.fun(:,:,:,tmpqi(1),tmpqi(2));
      if opt.hasmsk
        tmpmask = opt.msk(:,:,:,tmpqi(1),tmpqi(2));
      end
    end
    
    plotoptions = {'transform', eye(4), 'location', opt.ijk, ...
      'style', 'subplot', 'parents', tmph, 'update', opt.update, ...
      'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax]};
    if opt.hasmsk
      plotoptions = cat(2, plotoptions, {'datmask', tmpmask, 'opacitylim', [opt.opacmin opt.opacmax]});
    elseif opt.hasbackground
      % there's a background, in the absence of the mask, the fun should be
      % plotted with an opacity value of 0.5
      plotoptions = cat(2, plotoptions, {'datmask', ones(size(tmpfun))./2, 'opacitylim', [0 1]});
    end
    if opt.hasbackground
      % the background should always be added, independent of the mask
      plotoptions = cat(2, plotoptions, {'background', opt.background, 'maskstyle', 'colormix'});
    end
    ft_plot_ortho(tmpfun, plotoptions{:});
    % After the first call, the handles to the functional surfaces exist.
    % Create a variable containing these, and sort according to the parents.
    opt.funhandles = findobj(opt.handlesfigure, 'type', 'surface');
    opt.funtag     = get(opt.funhandles, 'tag');
    opt.funhandles = opt.funhandles(~strcmp('ana', opt.funtag));
    for i=1:length(opt.funhandles)
      opt.parenttag{i} = get(get(opt.funhandles(i), 'parent'), 'tag');
    end
    [i1,i2,i3] = intersect(opt.parenttag, {'ik' 'jk' 'ij'});
    opt.funhandles = opt.funhandles(i3(i2)); % seems like swapping the order
    opt.funhandles = opt.funhandles(:)';
    set(opt.funhandles, 'tag', 'fun');
    
    if ~opt.hasmsk && opt.hasfun && opt.hasana
      set(opt.funhandles(1), 'facealpha',0.5);
      set(opt.funhandles(2), 'facealpha',0.5);
      set(opt.funhandles(3), 'facealpha',0.5);
    end

  else
    if isequal(opt.funcolormap, 'rgb')
      tmpfun = opt.fun;
      if opt.hasmsk
        tmpmask = opt.msk;
      end
    else
      tmpqi  = [opt.qi 1];
      tmpfun = opt.fun(:,:,:,tmpqi(1),tmpqi(2));
      if opt.hasmsk
        tmpmask = opt.msk(:,:,:,tmpqi(1),tmpqi(2));
      end
    end
    
    plotoptions = {'transform', eye(4), 'location', opt.ijk, ...
        'style', 'subplot', 'surfhandle', opt.funhandles, 'update', opt.update, ...
        'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax]};
    if opt.hasmsk
      plotoptions = cat(2, plotoptions, {'datmask', tmpmask, 'opacitylim', [opt.opacmin opt.opacmax]});
    elseif opt.hasbackground
      % there's a background, in the absence of the mask, the fun should be
      % plotted with an opacity value of 0.5
      plotoptions = cat(2, plotoptions, {'datmask', ones(size(tmpfun))./2, 'opacitylim', [0 1]});
    end
    if opt.hasbackground
      % the background should always be added, independent of the mask
      plotoptions = cat(2, plotoptions, {'background', opt.background, 'maskstyle', 'colormix'});
    end
    ft_plot_ortho(tmpfun, plotoptions{:});
  end
end
set(opt.handlesaxes(1), 'Visible', opt.axis);
set(opt.handlesaxes(2), 'Visible', opt.axis);
set(opt.handlesaxes(3), 'Visible', opt.axis);

if opt.init
  if showcoordsys
    % add L/R label in the relevant panels
    ijk = 'ijk';
    [ind_ijk, ind_left] = find(strcmp(dirijk, 'left'));

    lr_tag = ijk(ind_ijk);
    lr_dir = [ind_left 3-ind_left];
    lr_str = 'LR';
    lr_str = lr_str(lr_dir);

    % by construction the handlesAxes 1/2/3 are ordered according to 'ik', 'jk', 'ij'
    for k = 1:3
      tag = get(opt.handlesaxes(k), 'Tag');
      if contains(tag, lr_tag)
        textcoord = [0 0 0; 0 0 0];
        [x, y]    = find(tag(:)==ijk);
        if find(tag==lr_tag)==2
          y = flip(y);
        end
        for kk = 1:2
          switch y(kk)
            case 1
              lims = get(opt.handlesaxes(k), 'Xlim');
            case 2
              lims = get(opt.handlesaxes(k), 'Ylim');
            case 3
              lims = get(opt.handlesaxes(k), 'Zlim');
          end
          if kk==1
            textcoord(1, y(kk)) = lims(1) + 14;
            textcoord(2, y(kk)) = lims(2) - 24;
          else
            textcoord(1, y(kk)) = lims(1) + 14;
            textcoord(2, y(kk)) = lims(1) + 14;
          end
        end
        last = setdiff(1:3, y);
        switch last
          case 1
            lims = get(opt.handlesaxes(k), 'Xlim');
          case 2
            lims = get(opt.handlesaxes(k), 'Ylim');
          case 3
            lims = get(opt.handlesaxes(k), 'Zlim');
        end
        if k==1
          lastval = -1;
        else
          lastval = lims(2)+1;
        end
        textcoord(1, last) = lastval;
        textcoord(2, last) = lastval;
        text(opt.handlesaxes(k), textcoord(:,1), textcoord(:,2), textcoord(:,3), ...
          lr_str(:), 'color', 'w', 'fontweight', 'bold');

      else
        continue;
      end
    end
  end
end


if opt.hasfreq && opt.hastime && opt.hasfun
  h4 = subplot(2,2,4);
  tmpdat = double(shiftdim(opt.fun(xi,yi,zi,:,:),3));
  % uimagesc is in external/fileexchange
  ft_hastoolbox('fileexchange', 1);
  uimagesc(double(functional.time), double(functional.freq), tmpdat); axis xy;
  xlabel('time'); ylabel('freq');
  set(h4, 'tag', 'TF1');
  caxis([opt.fcolmin opt.fcolmax]);
elseif opt.hasfreq && opt.hasfun
  h4 = subplot(2,2,4);
  plot(functional.freq, shiftdim(opt.fun(xi,yi,zi,:),3)); xlabel('freq');
  axis([functional.freq(1) functional.freq(end) opt.fcolmin opt.fcolmax]);
  set(h4, 'tag', 'TF2');
elseif opt.hastime && opt.hasfun
  h4 = subplot(2,2,4);
  plot(functional.time, shiftdim(opt.fun(xi,yi,zi,:),3)); xlabel('time');
  set(h4, 'tag', 'TF3', 'xlim',functional.time([1 end]), 'ylim',[opt.fcolmin opt.fcolmax], 'layer', 'top');
elseif strcmp(opt.colorbar,  'yes') && ~isfield(opt, 'hc')
  if opt.hasfun
    % vectorcolorbar = linspace(fscolmin, fcolmax,length(cfg.funcolormap));
    % imagesc(vectorcolorbar,1,vectorcolorbar);ft_colormap(cfg.funcolormap);
    % use a normal MATLAB colorbar, attach it to the invisible 4th subplot
    try
      caxis([opt.fcolmin opt.fcolmax]);
    end
    
    opt.hc = colorbar;
    set(opt.hc, 'location', 'southoutside');
    set(opt.hc, 'position',[0.06+0.06+opt.h1size(1) 0.06-0.06+opt.h3size(2) opt.h2size(1) 0.06]);
    ylabel(opt.hc, opt.colorbartext);
    try
      set(opt.hc, 'XLim', [opt.fcolmin opt.fcolmax]);
    end
  else
    ft_warning('no colorbar possible without functional data');
  end
end

if ~((opt.hasfreq && numel(functional.freq)>1) || opt.hastime)
  if opt.init
    ht = subplot('position',[0.06+0.06+opt.h1size(1) 0.06 opt.h2size(1) opt.h3size(2)]);
    set(ht, 'visible', 'off');
    opt.ht1=text(0,0.65,str1);
    opt.ht2=text(0,0.5,str2);
    opt.ht3=text(0,0.4,str4);
    opt.ht4=text(0,0.3,str3);
    opt.ht5=text(0,0.2,['atlas label: ' lab]);
  else
    set(opt.ht1, 'string',str1);
    set(opt.ht2, 'string',str2);
    set(opt.ht3, 'string',str4);
    set(opt.ht4, 'string',str3);
    set(opt.ht5, 'string',['atlas label: ' lab]);
  end
end

% make the last current axes current again
sel = findobj('type', 'axes', 'tag',tag);
if ~isempty(sel)
  set(opt.handlesfigure, 'currentaxes', sel(1));
end
if opt.crosshair
  if opt.init
    hch1 = ft_plot_crosshair([xi 1 zi], 'parent', opt.handlesaxes(1));
    hch3 = ft_plot_crosshair([xi yi opt.dim(3)], 'parent', opt.handlesaxes(3));
    hch2 = ft_plot_crosshair([opt.dim(1) yi zi], 'parent', opt.handlesaxes(2));
    opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  else
    ft_plot_crosshair([xi 1 zi], 'handle', opt.handlescross(1, :));
    ft_plot_crosshair([opt.dim(1) yi zi], 'handle', opt.handlescross(2, :));
    ft_plot_crosshair([xi yi opt.dim(3)], 'handle', opt.handlescross(3, :));
  end
end

if opt.init
  opt.init = false;
  setappdata(h, 'opt', opt);
end

set(h, 'currentaxes', curr_ax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');
tag     = get(curr_ax, 'tag');

if isempty(key)
  % this happens if you press the apple key
  key = '';
end

% the following code is largely shared by FT_SOURCEPLOT, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN, FT_MESHREALIGN, FT_ELECTRODEPLACEMENT
switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing
    
  case '1'
    subplot(opt.handlesaxes(1));
    
  case '2'
    subplot(opt.handlesaxes(2));
    
  case '3'
    subplot(opt.handlesaxes(3));
    
  case 'q'
    setappdata(h, 'opt', opt);
    cb_quit(h);
    
  case {'i' 'j' 'k' 'm' 28 29 30 31 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'} % TODO FIXME use leftarrow rightarrow uparrow downarrow
    % update the view to a new position
    if     strcmp(tag, 'ik') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    else
      % do nothing
    end
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case {43 'add' 'shift+equal'}  % + or numpad +
    % contrast scaling
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % reduce color scale range by 5%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)+cscalefactor;
    opt.clim(2) = opt.clim(2)-cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  case {45 'subtract' 'hyphen' 'shift+hyphen'} % - or numpad -
    % contrast scaling
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % increase color scale range by 5%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)-cscalefactor;
    opt.clim(2) = opt.clim(2)+cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
    
  otherwise
    % do nothing
    
end % switch key
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)

h   = getparent(h);
cb_getposition(h);

switch get(h, 'selectiontype')
  case 'normal'
    % just update to new position, nothing else to be done here
    cb_redraw(h);
  case 'alt'
    set(h, 'windowbuttonmotionfcn', @cb_tracemouse);
    cb_redraw(h);
  otherwise
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonrelease(h, eventdata)

set(h, 'windowbuttonmotionfcn', '');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_tracemouse(h, eventdata)

h   = getparent(h);
cb_getposition(h);
cb_redraw(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
pos     = mean(get(curr_ax, 'currentpoint'));

tag = get(curr_ax, 'tag');

if ~isempty(tag) && ~opt.init
  if strcmp(tag, 'ik')
    opt.ijk([1 3])  = round(pos([1 3]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'ij')
    opt.ijk([1 2])  = round(pos([1 2]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'jk')
    opt.ijk([2 3])  = round(pos([2 3]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF1')
    % timefreq
    opt.qi(2) = nearest(opt.functional.time, pos(1));
    opt.qi(1) = nearest(opt.functional.freq, pos(2));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF2')
    % freq only
    opt.qi  = nearest(opt.functional.freq, pos(1));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF3')
    % time only
    opt.qi  = nearest(opt.functional.time, pos(1));
    opt.update = [1 1 1];
  end
end
opt.ijk = min(opt.ijk(:)', opt.dim);
opt.ijk = max(opt.ijk(:)', [1 1 1]);

setappdata(h, 'opt', opt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)

% opt = getappdata(h, 'opt');
% opt.quit = true;
% setappdata(h, 'opt', opt);
% uiresume
delete(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

end
%%

function [omap] = alphamap_inourcode(param1, param2, param3)
global h
%ALPHAMAP - Set a figure's AlphaMap property
%
% ALPHAMAP(MATRIX)     - Set the current figure's AlphaMap property to MATRIX.
% ALPHAMAP('default')  - Set the AlphaMap to it's default value.
% ALPHAMAP('rampup')   - Create a linear alphamap_inourcode with increasing opacity.
% ALPHAMAP('rampdown') - Create a linear alphamap_inourcode with decreasing opacity.
% ALPHAMAP('vup')      - Create an alphamap_inourcode transparent in the center, and
%			 linearly increasing to the beginning and end.
% ALPHAMAP('vdown')    - Create an alphamap_inourcode opaque in the center, and
%			 linearly decreasing to the beginning and end.
% ALPHAMAP('increase') - Modify the alphamap_inourcode making it more opaque.
% ALPHAMAP('decrease') - Modify the alphamap_inourcode making it more transparent.
% ALPHAMAP('spin')     - Rotate the current alphamap_inourcode.
%
% ALPHAMAP(PARAM, LENGTH) - For Parameters which create new maps, create
%                        them with so they are LENGTH long.
% ALPHAMAP(CHANGE, DELTA) - For parameters which change the alphamap_inourcode, use
%                        DELTA as a parameter.
%
% ALPHAMAP(FIGURE,PARAM) - Set FIGURE's AlphaMap to some PARAMeter.
% ALPHAMAP(FIGURE,PARAM,LENGTH)
% ALPHAMAP(FIGURE,CHANGE)
% ALPHAMAP(FIGURE,CHANGE,DELTA)
%
% ALPHAMAP(AXES,PARAM) - If per axis alphamap_inourcode is supported, update axis AlphaMap to some PARAMeter,
% else set FIGURE that contains AXES AlphaMap.
% ALPHAMAP(AXES,PARAM,LENGTH)
% ALPHAMAP(AXES,CHANGE)
% ALPHAMAP(AXES,CHANGE,DELTA)
%
% AMAP=ALPHAMAP         - Fetch the current alphamap_inourcode
% AMAP=ALPHAMAP(FIGURE) - Fetch the current alphamap_inourcode from FIGURE.
% AMAP=ALPHAMAP(AXES) - If AXES alphamap_inourcode is supported, fetch AXES alphamap_inourcode else return alphamap_inourcode of FIGURE
% containing AXES.
% AMAP=ALPHAMAP(PARAM)  - Return the alphamap_inourcode based on PARAM
% 			  without setting the property.
%
% See also ALPHA, ALIM, COLORMAP.

% MAPSTRINGS=ALPHAMAP('strings') - Return a list of strings which generate
%                         alphamaps.

% Copyright 1984-2017 The MathWorks, Inc.

import matlab.graphics.internal.*;
set_alphamap = 0;
delta=0;

if nargin > 0

    hMapObject = getMapContainer_inourcode(param1);
    if isempty(hMapObject)
        len = size(get(h,'AlphaMap'),2);
    else
        len = size(get(hMapObject,'AlphaMap'),2);
    end
    
    if isscalar(param1) && ~isempty(hMapObject)        
        if nargin > 1
            param1 = param2;
            set_alphamap = 1;
            if nargin > 2
                len = param3;
                delta = param3;
            end
        else
            omap = get(hMapObject,'AlphaMap');
        end
    else
        hMapObject = h;
        if nargin > 0
            set_alphamap = 1;
            if nargin > 1
                len = param2;
                delta = param2;
            end
        else
            omap = get(hMapObject,'AlphaMap');
        end
    end
else
    hMapObject = h;
    omap = get(hMapObject,'AlphaMap');
    return
end

if isCharOrString(len)
    len = eval(len);
end
if isCharOrString(delta)
    delta = eval(delta);
end

if set_alphamap
    if isCharOrString(param1)
        switch param1
            case 'strings'
                map = { 'rampup' 'rampdown' 'vup' 'vdown' };
                set_alphamap = 0;
            case 'rampup'
                map = linspace(0, 1, len);
            case 'rampdown'
                map = linspace(1, 0, len);
            case 'vup'
                map = [linspace(0, 1, ceil(len/2)) linspace(1, 0, floor(len/2))];
            case 'vdown'
                map = [linspace(1, 0, ceil(len/2)) linspace(0, 1, floor(len/2))];
            case 'increase'
                map = get(hMapObject,'AlphaMap');
                if delta == 0
                    delta = .1;
                end
                map = map + delta;
                map(map > 1) = 1;
            case 'decrease'
                map = get(hMapObject,'AlphaMap');
                if delta == 0
                    delta = .1;
                end
                map = map - delta;
                map(map < 0) = 0;
            case 'spin'
                map = get(hMapObject,'AlphaMap');
                if delta == 0
                    delta = 1;
                end
                if delta > 0
                    map = [ map(delta+1:end) map(1:delta) ];
                elseif delta < 0
                    delta = - delta;
                    map = [ map(end-delta:end) map(1:end-delta-1) ];
                end
            case 'default'
                hFig = ancestor(hMapObject,'figure');
                map = get(hFig,'defaultfigureAlphamap');
            otherwise
                error(message('MATLAB:alphamap_inourcode:UnknownSpecifier'));
        end
    else
        map = param1;
    end
    
    if set_alphamap
        if nargout == 1
            omap = map;
        else
            set(hMapObject,'AlphaMap',map);
        end
    else
        omap = map;
    end
    
else
    omap = get(hMapObject,'AlphaMap');
end
end


%%

function mapContainer = getMapContainer_inourcode(obj)
global h
% getMapContainer(OBJECT), returns OBJECT
% If OBJECT has the AlphaMap/ColorMap property
% If OBJECT doesn't have AlphaMap/ColorMap as it's property but has a ColorSpace, then the ColorSpace is returned
% If OBJECT doesn't have either of the above,  EMPTY is returned.

% Does OBJ have a property called ColorMap or AlphaMap?
if isscalar(obj) && ( isprop(obj,'Colormap') || isprop(obj,'Alphamap') )
    mapContainer = obj;
    return;
end

% Does OBJ have a valid ColorSpace?
% ColorSpaces have the AlphaMap/ColorMap property,
if isscalar(obj) && isprop(obj, 'ColorSpace') && isa(get(obj,'ColorSpace'), 'matlab.graphics.axis.colorspace.ColorSpace')
    mapContainer = get(obj,'ColorSpace');
    return;
end

% object does not contain a map return empty;
mapContainer = [];
end
%%

function [hs] = ft_plot_mesh_inourcode(mesh, varargin)
global h
global ax
% FT_PLOT_MESH visualizes a surface or volumetric mesh, for example with the cortical
% folding of the brain, or the scalp surface of the head. Surface meshes are
% described by triangles and consist of a structure with the fields "pos" and "tri".
% Volumetric meshes are described with tetraheders or hexaheders and have the fields
% "pos" and "tet" or "hex".
%
% Use as
%   ft_plot_mesh(mesh, ...)
% or if you only want to plot the 3-D vertices
%   ft_plot_mesh(pos, ...)
%
% Optional arguments should come in key-value pairs and can include
%   'facecolor'    = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx3 or Nx1 array where N is the number of faces
%   'vertexcolor'  = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx3 or Nx1 array where N is the number of vertices
%   'edgecolor'    = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%   'faceindex'    = true or false
%   'vertexindex'  = true or false
%   'facealpha'    = transparency, between 0 and 1 (default = 1)
%   'edgealpha'    = transparency, between 0 and 1 (default = 1)
%   'surfaceonly'  = true or false, plot only the outer surface of a hexahedral or tetrahedral mesh (default = false)
%   'vertexmarker' = character, e.g. '.', 'o' or 'x' (default = '.')
%   'vertexsize'   = scalar or vector with the size for each vertex (default = 10)
%   'unit'         = string, convert to the specified geometrical units (default = [])
%   'axes'          = boolean, whether to plot the axes of the 3D coordinate system (default = false)
%   'maskstyle',   = 'opacity' or 'colormix', if the latter is specified, opacity masked color values
%                    are converted (in combination with a background color) to RGB. This bypasses
%                    openGL functionality, which behaves unpredictably on some platforms (e.g. when
%                    using software opengl)
%
% If you don't want the faces, edges or vertices to be plotted, you should specify the color as 'none'.
%
% Example
%   [pos, tri] = mesh_sphere(162);
%   mesh.pos = pos;
%   mesh.tri = tri;
%   ft_plot_mesh(mesh, 'facecolor', 'skin', 'edgecolor', 'none')
%   camlight
%
% You can plot an additional contour around specified areas using
%   'contour'           = inside of contour per vertex, either 0 or 1
%   'contourcolor'      = string, color specification
%   'contourlinestyle'  = string, line specification
%   'contourlinewidth'  = number
%
% See also FT_PLOT_HEADSHAPE, FT_PLOT_HEADMODEL, TRIMESH, PATCH

% Copyright (C) 2009, Cristiano Micheli
% Copyright (C) 2009-2022, Robert Oostenveld
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

% rename pnt into pos
mesh = fixpos(mesh);

if ~isstruct(mesh) && isnumeric(mesh) && size(mesh,2)==3
  % the input seems like a list of points, convert into something that resembles a mesh
  mesh = struct('pos', mesh);
end

% the input is a structure, but might also be a struct-array
if numel(mesh)>1
  % plot each of the boundaries
  for i=1:numel(mesh)
    ft_plot_mesh_inourcode(mesh(i), varargin{:})
  end
  return
end

% get the optional input arguments
vertexcolor  = ft_getopt(varargin, 'vertexcolor');
if isfield(mesh, 'tri') && size(mesh.tri,1)>10000
  facecolor    = ft_getopt(varargin, 'facecolor',   'cortex_light');
  edgecolor    = ft_getopt(varargin, 'edgecolor',   'none');
else
  facecolor    = ft_getopt(varargin, 'facecolor',   'white');
  edgecolor    = ft_getopt(varargin, 'edgecolor',   'k');
end
% faceindex=0;
% vertexindex
% vertexsize
% vertexmarker
% facealpha
% edgealpha
% edgelinewidth
% material_
% tag
% surfaceonly
% unit
% axes_
% clim
% alphalim
% alphamapping
% cmap
% maskstyle, contour, contourcolor, contourlinestyle, contourlinewidth
% 
% 

















faceindex     = ft_getopt(varargin, 'faceindex',   false);
vertexindex   = ft_getopt(varargin, 'vertexindex', false);
vertexsize    = ft_getopt(varargin, 'vertexsize',  10);
vertexmarker  = ft_getopt(varargin, 'vertexmarker', '.');
facealpha     = ft_getopt(varargin, 'facealpha',   1);
edgealpha     = ft_getopt(varargin, 'edgealpha',   1);
edgelinewidth = ft_getopt(varargin, 'edgelinewidth', .5);
material_     = ft_getopt(varargin, 'material');        % note the underscore, there is also a material function
tag           = ft_getopt(varargin, 'tag',         '');
surfaceonly   = ft_getopt(varargin, 'surfaceonly');     % default is handled below
unit          = ft_getopt(varargin, 'unit');
axes_         = ft_getopt(varargin, 'axes', false);     % do not confuse with built-in function
clim          = ft_getopt(varargin, 'clim');
alphalim      = ft_getopt(varargin, 'alphalim');
alphamapping  = ft_getopt(varargin, 'alphamap_inourcode', 'rampup');
cmap          = ft_getopt(varargin, 'colormap');
maskstyle     = ft_getopt(varargin, 'maskstyle', 'opacity');
contour       = ft_getopt(varargin, 'contour',   []);

contourcolor      = ft_getopt(varargin, 'contourcolor',     'k');
contourlinewidth  = ft_getopt(varargin, 'contourlinewidth', 3);
contourlinestyle  = ft_getopt(varargin, 'contourlinestyle', '-');

haspos   = isfield(mesh, 'pos');   % vertices
hastri   = isfield(mesh, 'tri');   % triangles   as a Mx3 matrix with vertex indices
hastet   = isfield(mesh, 'tet');   % tetraheders as a Mx4 matrix with vertex indices
hashex   = isfield(mesh, 'hex');   % hexaheders  as a Mx8 matrix with vertex indices
hasline  = isfield(mesh, 'line');  % line segments in 3-D
haspoly  = isfield(mesh, 'poly');  % polygons describing a surface in 3-D
hascolor = isfield(mesh, 'color'); % color code for vertices

if isempty(surfaceonly)
  if hastet
    ft_warning('only visualizing the outer surface of the tetrahedral mesh, see the "surfaceonly" option')
    surfaceonly = true;
  elseif hashex
    ft_warning('only visualizing the outer surface of the hexahedral mesh, see the "surfaceonly" option')
    surfaceonly = true;
  else
    surfaceonly = false;
  end
end

if ~isempty(unit)
  mesh = ft_convert_units(mesh, unit);
end

if surfaceonly
  mesh = mesh2edge(mesh);
  % update the flags that indicate which surface/volume elements are present
  hastri   = isfield(mesh, 'tri');  % triangles   as a Mx3 matrix with vertex indices
  hastet   = isfield(mesh, 'tet');  % tetraheders as a Mx4 matrix with vertex indices
  hashex   = isfield(mesh, 'hex');  % hexaheders  as a Mx8 matrix with vertex indices
  haspoly  = isfield(mesh, 'poly'); % polygons
end

% convert string into boolean values
faceindex   = istrue(faceindex);   % yes=view the face number
vertexindex = istrue(vertexindex); % yes=view the vertex number

if isempty(vertexcolor)
  if haspos && hascolor && (hastri || hastet || hashex || hasline || haspoly)
    vertexcolor = mesh.color;
  elseif haspos && (hastri || hastet || hashex || hasline || haspoly)
    vertexcolor ='none';
  else
    vertexcolor ='k';
  end
end

% there are various ways of specifying that this should not be plotted
if isequal(vertexcolor, 'false') || isequal(vertexcolor, 'no') || isequal(vertexcolor, 'off') || isequal(vertexcolor, false)
  vertexcolor = 'none';
end
if isequal(facecolor, 'false') || isequal(facecolor, 'no') || isequal(facecolor, 'off') || isequal(facecolor, false)
  facecolor = 'none';
end
if isequal(edgecolor, 'false') || isequal(edgecolor, 'no') || isequal(edgecolor, 'off') || isequal(edgecolor, false)
  edgecolor = 'none';
end

% color management
if ischar(vertexcolor) && exist([vertexcolor '.m'], 'file')
  vertexcolor = feval(vertexcolor);
elseif ischar(vertexcolor) && ismember(vertexcolor, htmlcolors)
  vertexcolor = htmlcolors(vertexcolor);
elseif ischar(vertexcolor) && isequal(vertexcolor, 'curv') % default of ft_sourceplot method surface
  if isfield(mesh, 'curv')
    cortex_light = feval('cortex_light');
    cortex_dark  = feval('cortex_dark');
    % the curvature determines the color of gyri and sulci
    vertexcolor = mesh.curv(:) * cortex_dark + (1-mesh.curv(:)) * cortex_light;
  else
    cortex_light = [199 194 169]/255;
    vertexcolor = repmat(cortex_light, size(mesh.pos,1), 1);
    ft_warning('no curv field present in the mesh structure, using cortex_light as vertexcolor')
  end
end

if ischar(facecolor) && exist([facecolor '.m'], 'file')
  facecolor = [1,1,1];
elseif ischar(facecolor) && ismember(facecolor, htmlcolors)
  facecolor = htmlcolors(facecolor);
end

if ischar(edgecolor) && exist([edgecolor '.m'], 'file')
  edgecolor = 'none';
elseif ischar(edgecolor) && ismember(edgecolor, htmlcolors)
  edgecolor = htmlcolors(edgecolor);
end

% everything is added to the current figure
% ax.NextPlot = 'add';
% holdflag = ishold;
% if ~holdflag
%   hold on
% end

if isfield(mesh, 'pos')
  % this is assumed to reflect 3-D vertices
  pos = mesh.pos;
elseif isfield(mesh, 'prj')
  % this happens sometimes if the 3-D vertices are projected to a 2-D plane
  pos = mesh.prj;
else
  ft_error('no vertices found');
end

if isempty(pos)
  hs=[];
  return
end

if hastri+hastet+hashex+hasline+haspoly>1
  ft_error('cannot deal with simultaneous triangles, tetraheders and/or hexaheders')
end

if hastri
  tri = mesh.tri;
elseif haspoly
  % these are treated just like triangles
  tri = mesh.poly;
elseif hastet
  % represent the tetraeders as the four triangles
  tri = [
    mesh.tet(:,[1 2 3]);
    mesh.tet(:,[2 3 4]);
    mesh.tet(:,[3 4 1]);
    mesh.tet(:,[4 1 2])];
  % or according to SimBio:  (1 2 3), (2 4 3), (4 1 3), (1 4 2)
  % there are shared triangles between neighbouring tetraeders, remove these
  tri = unique(tri, 'rows');
elseif hashex
  % represent the hexaheders as a collection of 6 patches
  tri = [
    mesh.hex(:,[1 2 3 4]);
    mesh.hex(:,[5 6 7 8]);
    mesh.hex(:,[1 2 6 5]);
    mesh.hex(:,[2 3 7 6]);
    mesh.hex(:,[3 4 8 7]);
    mesh.hex(:,[4 1 5 8]);
    ];
  % there are shared faces between neighbouring hexaheders, remove these
  tri = unique(tri, 'rows');
else
  tri = [];
end

if hasline
  line = mesh.line;
else
  line = [];
end

if haspos
  if ~isempty(tri)
    hs = patch(ax, 'Vertices', pos, 'Faces', tri);
  elseif ~isempty(line)
    hs = patch(ax, 'Vertices', pos, 'Faces', line);
  else
    hs = patch(ax, 'Vertices', pos, 'Faces', []);
  end
  %set(hs, 'FaceColor', facecolor);
  set(hs, 'EdgeColor', edgecolor);
  set(hs, 'LineWidth', edgelinewidth);
  set(hs, 'tag', tag);
end

if ~isempty(material_)
  material(material_); % dull, shiny or default
end

% the vertexcolor can be specified either as a RGB color for each vertex, or as a single value at each vertex
% the facecolor can be specified either as a RGB color for each triangle, or as a single value at each triangle
% if there are triangles, the vertexcolor is used for linear interpolation over the patches
vertexpotential = ~isempty(tri) && ~ischar(vertexcolor) && (size(pos,1)==numel(vertexcolor) || size(pos,1)==size(vertexcolor,1) && (size(vertexcolor,2)==1 || size(vertexcolor,2)==3));
facepotential   = ~isempty(tri) && ~ischar(facecolor  ) && (size(tri,1)==numel(facecolor  ) || size(tri,1)==size(facecolor  ,1) && (size(facecolor  ,2)==1 || size(facecolor,  2)==3));

switch maskstyle
  case 'opacity'
    % if both vertexcolor and facecolor are numeric arrays, let the vertexcolor prevail
    if vertexpotential
      % vertexcolor is an array with number of elements equal to the number of vertices
      set(hs, 'FaceVertexCData', vertexcolor, 'FaceColor', 'interp');
      if numel(vertexcolor)==size(pos,1)
        if ~isempty(clim), set(ax, 'clim', clim); end
        if ~isempty(cmap), ft_colormap_inourcode(cmap); end
      end
    elseif facepotential
      set(hs, 'FaceVertexCData', facecolor, 'FaceColor', 'flat');
      if numel(facecolor)==size(tri,1)
        if ~isempty(clim), set(ax, 'clim', clim); end
        if ~isempty(cmap), ft_colormap(cmap); end
      end
    else
      % the color is indicated as a single character or as a single RGB triplet
      set(hs, 'FaceColor', facecolor);
    end
    
    % facealpha is a scalar, or an vector matching the number of vertices
    if size(pos,1)==numel(facealpha)
      set(hs, 'FaceVertexAlphaData', facealpha);
      set(hs, 'FaceAlpha', 'interp');
    elseif ~isempty(pos) && numel(facealpha)==1 && facealpha~=1
      % the default is 1, so that does not have to be set
      set(hs, 'FaceAlpha', facealpha);
    end
    
    if edgealpha~=1
      % the default is 1, so that does not have to be set
      set(hs, 'EdgeAlpha', edgealpha);
    end
    
    if ~(all(facealpha==1) && edgealpha==1)
      if ~isempty(alphalim)
        alim(ax, alphalim);
      end
      alphamap_inourcode(alphamapping);
    end
    
  case 'colormix'
    % ensure facecolor to be 1x3
    assert(isequal(size(facecolor),[1 3]), 'facecolor should be 1x3');
    
    % ensure facealpha to be nvertex x 1
    if numel(facealpha)==1
      facealpha = repmat(facealpha, size(pos,1), 1);
    end
    assert(isequal(numel(facealpha),size(pos,1)), 'facealpha should be %dx1', size(pos,1));
    
    bgcolor = repmat(facecolor, [numel(vertexcolor) 1]);
    rgb     = bg_rgba2rgb(bgcolor, vertexcolor, cmap, clim, facealpha, alphamapping, alphalim);
    set(hs, 'FaceVertexCData', rgb, 'facecolor', 'interp');
    if ~isempty(clim); caxis(clim); end % set colorbar scale to match [fcolmin fcolmax]
end

if ~isempty(contour)
  if ~iscell(contour), contour = {contour}; end
  if ~iscell(contourlinestyle), contourlinestyle = {contourlinestyle}; end
  
  if ischar(contourcolor)
    if numel(contour)>numel(contourcolor)
      contourcolor = repmat(contourcolor(:), [numel(contour) 1]);
    else
      contourcolor = contourcolor(:);
    end
  end
  if size(contourcolor,2)==3 && numel(contour)>size(contourcolor,1), contourcolor = repmat(contourcolor, [numel(contour) 1] ); end
  if numel(contour)>numel(contourlinewidth), contourlinewidth = repmat(contourlinewidth, [1 numel(contour)]); end
  if numel(contour)>numel(contourlinestyle), contourlinestyle = repmat(contourlinestyle, [1 numel(contour)]); end
  
  for m = 1:numel(contour)
    C    = full(triangle2connectivity(tri));
    clus = findcluster(contour{m},C,0);
    
    for cl = 1:max(clus)
      idxcl = find(clus==cl);
      [xbnd, ybnd, zbnd] = extract_contour(pos,tri,idxcl,contour{m});
      
      % draw each individual line segment of the intersection
      p = [];
      for i = 1:length(xbnd)
        p(i) = patch(xbnd(i,:)', ybnd(i,:)', zbnd(i,:)',NaN);
      end
      set(p(:), 'EdgeColor', contourcolor(m,:), 'LineWidth', contourlinewidth(m), 'LineStyle', contourlinestyle{m});
    end
  end
end

if faceindex
  % plot the triangle indices (numbers) at each face
  for face_indx=1:size(tri,1)
    str = sprintf('%d', face_indx);
    tri_x = (pos(tri(face_indx,1), 1) +  pos(tri(face_indx,2), 1) +  pos(tri(face_indx,3), 1))/3;
    tri_y = (pos(tri(face_indx,1), 2) +  pos(tri(face_indx,2), 2) +  pos(tri(face_indx,3), 2))/3;
    tri_z = (pos(tri(face_indx,1), 3) +  pos(tri(face_indx,2), 3) +  pos(tri(face_indx,3), 3))/3;
    h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hs  = [hs; h];
  end
end

if ~isequal(vertexcolor, 'none') && ~vertexpotential
  % plot the vertices as points
  
  if isempty(vertexcolor)
    % use black for all points
    if isscalar(vertexsize)
      if size(pos,2)==2
        hs = plot(pos(:,1), pos(:,2), ['k' vertexmarker]);
      else
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), ['k' vertexmarker]);
      end
      set(hs, 'MarkerSize', vertexsize);
    else
      if size(pos,2)==2
        for i=1:size(pos,1)
          hs = plot(pos(i,1), pos(i,2), ['k' vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      else
        for i=1:size(pos,1)
          hs = plot3(pos(i,1), pos(i,2), pos(i,3), ['k' vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==1
    % one color for all points
    if isscalar(vertexsize)
      if size(pos,2)==2
        hs = plot(pos(:,1), pos(:,2), [vertexcolor vertexmarker]);
      else
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), [vertexcolor vertexmarker]);
      end
      set(hs, 'MarkerSize', vertexsize);
    else
      if size(pos,2)==2
        for i=1:size(pos,1)
          hs = plot(pos(i,1), pos(i,2), [vertexcolor vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      else
        for i=1:size(pos,1)
          hs = plot3(pos(i,1), pos(i,2), pos(i,3), [vertexcolor vertexmarker]);
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ischar(vertexcolor) && numel(vertexcolor)==size(pos,1)
    % one color for each point
    if size(pos,2)==2
      for i=1:size(pos,1)
        hs = plot(pos(i,1), pos(i,2), [vertexcolor(i) vertexmarker]);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize);
        else
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    else
      for i=1:size(pos,1)
        hs = plot3(pos(i,1), pos(i,2), pos(i,3), [vertexcolor(i) vertexmarker]);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize);
        else
          set(hs, 'MarkerSize', vertexsize(i));
        end
      end
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==1
    % one RGB color for all points
    if size(pos,2)==2
      hs = plot(pos(:,1), pos(:,2), vertexmarker);
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    else
      hs = plot3(pos(:,1), pos(:,2), pos(:,3), vertexmarker);
      set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor);
    end
    
  elseif ~ischar(vertexcolor) && size(vertexcolor,1)==size(pos,1) && size(vertexcolor,2)==3
    % one RGB color for each point
    if size(pos,2)==2
      for i=1:size(pos,1)
        hs = plot(pos(i,1), pos(i,2), vertexmarker);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
        else
          set(hs, 'MarkerSize', vertexsize(i), 'MarkerEdgeColor', vertexcolor(i,:));
        end
      end
    else
      for i=1:size(pos,1)
        hs = plot3(pos(i,1), pos(i,2), pos(i,3), vertexmarker);
        if isscalar(vertexsize)
          set(hs, 'MarkerSize', vertexsize, 'MarkerEdgeColor', vertexcolor(i,:));
        else
          set(hs, 'MarkerSize', vertexsize(i), 'MarkerEdgeColor', vertexcolor(i,:));
        end
      end
    end
    
  else
    ft_error('Unknown color specification for the vertices');
  end
  
end % plotting the vertices as points

if vertexindex
  % plot the vertex indices (numbers) at each node
  for node_indx=1:size(pos,1)
    str = sprintf('%d', node_indx);
    if size(pos, 2)==2
      h = text(pos(node_indx, 1), pos(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
      h = text(pos(node_indx, 1), pos(node_indx, 2), pos(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    hs = [hs; h];
  end
end

set(ax,'Visible','off');
set(get(ax,'Title'),'Visible','on');
set(ax,'CameraViewAngle',get(ax,'CameraViewAngle'));
set(ax,'PlotBoxAspectRatio',get(ax,'PlotBoxAspectRatio'));
set(ax,'DataAspectRatio',get(ax,'DataAspectRatio'));
% LocSetEqual(ax,pbarlimit);           

% axis off
% axis vis3d
% axis equal

if istrue(axes_)
  % plot the 3D axes, this depends on the units and coordsys
  ft_plot_axes(mesh);
end

if isfield(mesh, 'coordsys')
  % add a context sensitive menu to change the 3d viewpoint to top|bottom|left|right|front|back
  menu_viewpoint(h, mesh.coordsys)
end

% if ~holdflag
%   hold off
% end
global counter
if (counter ==2)
    
    drawnow;
end

if ~nargout
  clear hs
end
end
%%
%   skull             = [140  85  85]/255
%   cortex            = [255 213 119]/255;
%   cortex_light      = [199 194 169]/255;
%   cortex_dark       = [100  97  85]/255;
%   skin              = [249 223 192]/255;
%   skin_light        = [249 223 192]/255;
%   skin_medium_light = [225 194 158]/255;
%   skin_medium       = [188 142 106]/255;
%   skin_medium_dark  = [155 102	65]/255;
%   skin_dark         = [ 91  71  61]/255;






