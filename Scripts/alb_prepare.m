% ALB_PREPARE  Original daily GeoTIFFs -> analysis-ready yearly layers (GeoTIFF + CSV).
%
% Part 1 of the Alaska spring-albedo study runs as two scripts:
%   1. alb_prepare.m  (this file)  reads the ORIGINAL yearly GeoTIFFs, does all
%                                  preprocessing and writes one GeoTIFF per quantity
%                                  (one band per year) plus CSV logs.
%   2. alb_analysis.m              reads only those GeoTIFFs / CSVs and answers RQ1-RQ3.
% No .mat file is read or written anywhere.
%
% RESOLUTION (checked on the raw 2004, 2013 and 2023 files of the v1 archive,
% Sep 2026; the v3 albedo has more valid 1-km values - recheck these numbers)
%   Land cells with a valid spring transition in a year:
%       1 km 19-29 %    2 km 35-44 %    5 km 36-45 %    10 km 35-45 %
%   5-km SCD50 vs the median of its own 1-km pixels: median difference 0.0 d
%   (IQR -1.5 to +1.5 d). The spread of 1-km SCD50 inside a 5-km cell (SD 2-4 d,
%   partly retrieval noise) is ~20 % of the spatial variance. Synchrony length
%   scales are 100-300 km and ERA5-Land is ~9 km. Hence: 5 km is the primary grid;
%   10 km is produced in the same pass as a resolution check (C.aggs).
%
% INPUTS (under C.root; every raster is on the 1-km EPSG:3338 grid whose
% upper-left corner is (-706000, 2375000); the covariate files are one row and one
% column larger than the albedo files and the top-left crop aligns)
%   BlueSky_Daily_TIFFs_v4\BlueSky_albedo_daily_YYYY.tif   (GEE export v3, Oct 2026)
%       MCD43A3-based blue-sky albedo, int16, x 1e-4, one band per day,
%       raw 1..10000 valid (0 = fill, -32768 = nodata). Each 1-km value is the
%       area-weighted mean of the native 463-m pixels with BRDF quality <= 1
%       (full AND magnitude inversions), kept where they cover >= 50 % of the
%       pixel. Raw values above 10000 (albedo > 1) are counted in
%       prepare_log.csv and treated as missing.
%   BlueSky_Daily_QA_TIFFs_v4\BlueSky_pctfull_daily_YYYY.tif   (if C.useAlbQA)
%       companion of each albedo file: uint8, % of the value's area that comes
%       from FULL inversions (0..100; 255 = nodata); same grid and
%       NoData pattern; it must hold every window band of the albedo file
%       (matched by date; dates outside the window may differ, e.g. a partial
%       year whose two files were exported on different days). Gives mag_frac
%       and the full-inversion sensitivity run.
%   The v1 archive (BlueSky_Daily_TIFFs_Clipped: GEE default resampling, full
%   inversions only) still runs: point VAR.alb at it and set C.useAlbQA = false.
%   SnowCover_Daily_TIFFs\NDSI_SnowCover_MOD10A1_daily_YYYY.tif
%       MOD10A1 NDSI snow cover, int16, 0..100 valid (-> fraction), validation only
%   ERA5_Daily_TIFFs\ERA5_ssrd_daily_YYYY.tif        (only if C.useERA5)
%       ERA5-Land surface solar radiation downwards, daily sum, int16,
%       raw x 0.002 = MJ m-2 d-1 (export snippet at the end of this header)
%   topography\ALOS_AW3D30_topography.tif
%       float32 bands elevation_m, slope_deg, aspect_deg; band 1 finite = the
%       Alaska AOI used for every GEE export
%   ecoregions_simplified\ecoregions_simplified.shp   9 Level-2 ecoregions
%   Band dates are read from the GDAL band descriptions 'dYYYYMMDD' (tag 42112,
%   read from the file itself: imfinfo does not return it for BigTIFF, the GEE
%   export format); a file without them stops the run. Window days without a
%   band (MODIS outages, e.g. albedo 2001 days 175-176; NDSI 2001 days 167-183
%   and 2002 days 79-86 and 105) enter as days without observations; a note
%   lists them and prepare_log.csv counts them.
%
% PREPROCESSING, for every year and every grid (agg = 5 and 10 km)
%   1. Window: the bands whose date falls on days C.window of a 365-day calendar
%      (60 = 1 Mar ... 200 = 19 Jul in every year; in leap years 29 Feb is
%      skipped), so a day number means the same calendar date in all years.
%   2. Quality: a 1-km value is valid if its raw value lies in the variable's
%      valid range (VAR below). Raw values that are neither valid nor nodata/fill
%      are counted and written to prepare_log.csv. With C.useAlbQA, an albedo
%      value is also dropped when less than C.minPctFull % of its area comes
%      from full inversions (0 = keep all, the primary run; 100 = full
%      inversions only, the sensitivity run).
%   3. Aggregation: non-overlapping agg x agg blocks from the top-left corner of
%      the 1-km grid. Cell-day value = mean of the valid 1-km pixels; kept when at
%      least ceil(C.minFrac x agg^2) of them are valid (5 of 25 at 5 km).
%   4. Spring transition of each cell-year (numerics identical to alb_metrics.m):
%      Gaussian kernel smoother on the observed days (Nadaraya-Watson, bandwidth
%      6 d, smoothed value kept where the summed kernel weight is >= 0.8);
%      a_snow = 85th percentile of the smoothed curve on days 60-100,
%      a_floor = 15th percentile on days 165-200, amp = a_snow - a_floor;
%      scd75 / scd50 / scd25 = first downward crossing of a_floor + 0.75 / 0.50 /
%      0.25 amp (linear interpolation between smoothed days); dur = scd25 - scd75;
%      k_rel = 0.5 / dur (share of the amplitude lost per day in the central half).
%      QC: amp > 0.05, 0 < dur < 80 d and >= 20 valid days. qc code: 0 valid,
%      1 < 20 valid days, 2 a plateau window has no smoothed value, 3 amp <= 0.05
%      (no seasonal snow signal), 4 no crossing or duration outside (0, 80) d.
%      The same metric on NDSI snow cover gives ndsi_scd50 (validation).
%   5. Spring shortwave over days C.window, for cell-years that pass QC. Daily
%      albedo = the smoothed curve; days without a smoothed value are filled
%      (inside the record: linear in time; before the first smoothed day: a_snow;
%      after the last: a_floor). fill_frac = share of the incoming shortwave that
%      falls on filled days; the analysis decides which fill_frac is acceptable.
%          sw_in  = sum_d SWdown_d                  [MJ m-2]
%          sw_abs = sum_d SWdown_d (1 - albedo_d)   [MJ m-2]
%      SWdown = top-of-atmosphere insolation (always; Spencer 1971 / Iqbal 1983,
%      S0 = C.S0) and, if C.useERA5, the ERA5-Land all-sky surface flux (a cell-day
%      of SSRD = mean of the valid 1-km pixels, >= 1 pixel).
%   6. Static layers: ecoregion (1-km rasterisation of the shapefile, then the
%      code held by >= C.majorityFrac of the cell's 1-km pixels, else 0), land
%      fraction (share of 1-km pixels inside the AOI; a cell is analysed when
%      >= C.landFrac), latitude / longitude of the cell centre, mean and standard
%      deviation of elevation over the cell's AOI pixels; mean slope over the AOI
%      pixels where the DEM slope is defined (it is NaN on ~7 % of AOI pixels along
%      coasts and the AOI edge).
%
% ECOREGION GEOMETRY (verified with alb_ecoregions.m, Sep 2026)
%   The simplified shapefile carries NAD27 / Alaska Albers coordinates (Clarke
%   1866) under a NAD83 label. They are converted NAD27 Albers -> lat/lon (Clarke
%   1866) -> geocentric shift NAD27 -> NAD83 (dX -5, dY +135, dZ +172 m; NIMA
%   TR8350.2, Alaska) -> EPSG:3338. The bounding box is checked against the
%   original NAD27 extent. Rasterisation: pixel-centre rule, even-odd for holes.
%   Nine named places must fall in their expected region (E.checkPts).
%
% OUTPUTS  <C.outRoot>\grid_5km\  and  <C.outRoot>\grid_10km\
%   <quantity>.tif   one band per year; band k = year k of years.csv
%       scd75, scd50, scd25 [day, 365-day calendar], dur [d], k_rel [1/d],
%       a_snow, a_floor, amp [albedo], nobs [valid days], qc [code 0-4],
%       ndsi_scd50 [day], ndsi_qc [code 0-4],
%       sw_in_toa, sw_abs_toa [MJ m-2], fill_frac_toa [0-1],
%       sw_in_era5, sw_abs_era5 [MJ m-2], fill_frac_era5 [0-1]   (C.useERA5 only)
%       mag_frac [0-1]   (C.useAlbQA only) share of the observed area behind the
%           cell's albedo series that comes from magnitude inversions: over the
%           window cell-days kept by step 3, sum of (1 - pct_full/100) of the
%           valid 1-km values / their number. Written for every land cell with
%           a kept cell-day (not only qc = 0); 0 when C.minPctFull = 100.
%     NaN outside the land cells and, for the metrics, where qc > 0.
%   static_<name>.tif   ecoregion, land_frac, lat_deg, lon_deg, elev_mean_m,
%                       elev_sd_m, slope_mean_deg
%   years.csv (band <-> year), ecoregion_codes.csv, prepare_log.csv (one row
%   per year and grid; incl. alb_ / ndsi_window_days_missing and, with
%   C.useAlbQA, alb_mag_area_pct, alb_qa_dropped, alb_qa_mismatch),
%   settings.csv (every setting of this run)
%   <C.outRoot>\prepare_console.txt (everything printed during the run)
%   A folder that already holds layers made from other albedo files is not
%   overwritten (rename it, e.g. Prepared_v1, or change C.outRoot).
%
% FULL-INVERSION SENSITIVITY RUN (read by alb_analysis.m as variant 'fullinv')
%   Same settings, plus C.minPctFull = 100 and
%   C.outRoot = fullfile(C.root, 'Prepared_fullinv').
%
% RUN TIME  roughly 2-4 min per year (two files of 365 bands each; about 1-1.5 h for
%           2002-2025), about 50 % more with ERA5, somewhat more with the QA
%           files (a second, smaller file per year). Memory: ~2 GB.
% Requires the Mapping Toolbox (shaperead, readgeoraster, georasterinfo,
% geotiffinfo, geotiffwrite, maprefcells) and the Tiff class (base MATLAB).
%
% ERA5-LAND SSRD EXPORT (Google Earth Engine Code Editor) - reuse the export
% settings of your existing ERA5_t2m files; only the band, the scale and the
% file name change. The script checks grid, band dates and value range.
%   var aoi = /* the AOI geometry used for the other exports */;
%   for (var y = 2002; y <= 2025; y++) {
%     var col = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
%         .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y + 1, 1, 1))
%         .select('surface_solar_radiation_downwards_sum');           // J m-2 per UTC day
%     var names = col.aggregate_array('system:index')
%         .map(function (s) { return ee.String('d').cat(s); });         // 'dYYYYMMDD'
%     var img = col.map(function (im) { return im.divide(2000).round(); })  // raw x 0.002 = MJ m-2
%         .toBands().rename(names).toInt16().unmask(-32768).clip(aoi);
%     Export.image.toDrive({image: img, description: 'ERA5_ssrd_daily_' + y,
%         fileNamePrefix: 'ERA5_ssrd_daily_' + y, crs: 'EPSG:3338',
%         crsTransform: [1000, 0, -706000, 0, -1000, 2375000], dimensions: '2199x1834',
%         maxPixels: 1e10, formatOptions: {noData: -32768}});
%   }
% Note: ERA5 daily sums are UTC days (Alaska local midnight = 09 UTC); over a
% 141-day window this shifts the total by well under 1 %.

clear; clc;

%% ======================= USER SETTINGS =======================
C.root      = 'C:\Users\aymane.ahajjam\Desktop\Albedo';
C.outRoot   = fullfile(C.root, 'Prepared');   % full-inversion run: fullfile(C.root, 'Prepared_fullinv')
C.years     = 2000:2026;       % 2000-2001 are Terra-only; Aqua joined MCD43 in mid-2002
C.aggs      = [1];        %5 10   % km per grid cell; the first is the primary grid
C.window    = [60 200];        % days (365-day calendar): 1 Mar - 19 Jul
C.minFrac   = 0.2;             % cell-day kept if >= ceil(minFrac * agg^2) valid 1-km pixels
C.landFrac  = 0.5;             % cell analysed if >= this share of its 1-km pixels is in the AOI
C.majorityFrac = 0.5;          % ecoregion of a cell: code held by >= this share of its pixels
C.S0        = 1361;            % W m-2, total solar irradiance (Kopp & Lean 2011)
C.useERA5   = true;           % true once the ERA5_ssrd_daily_YYYY.tif files are in place
C.bandDates = 'tag';           % 'tag': dates from the GDAL band descriptions dYYYYMMDD (checked);
                               % 'dayOfYear': band k = day k of the year (use only if MATLAB cannot
                               % see the descriptions; the band count must then equal the year length)
C.useAlbQA   = true;           % read the QA companion files (VAR.albqa): mag_frac layer, C.minPctFull.
                               % false for albedo files without them (v1 archive)
C.minPctFull = 0;              % keep a 1-km albedo value only if >= this % of its area is full BRDF
                               % inversions: 0 = all (primary run), 100 = full inversions only
                               % (sensitivity run; set C.outRoot to Prepared_fullinv)

% spring-transition metric (alb_metrics.m values; all fields required)
C.opt = struct('bw', 6, 'minwt', 0.8, 'pre', [60 100], 'post', [165 200], ...
               'qhi', 85, 'qlo', 15, 'minamp', 0.05, 'minobs', 20);

% variables: folder, file pattern, scale to physical units, valid raw range
VAR.alb  = struct('dir', fullfile(C.root, 'BlueSky_Daily_TIFFs_v4'), ...
                  'pat', 'BlueSky_albedo_daily_%d.tif', 'scale', 1e-4, 'rawLo', 1, 'rawHi', 10000, ...
                  'fill', [0 -32768]);                    % albedo raw 0 is fill, not a black surface
VAR.albqa = struct('dir', fullfile(C.root, 'BlueSky_Daily_QA_TIFFs_v4'), ...
                  'pat', 'BlueSky_pctfull_daily_%d.tif', 'scale', 1, 'rawLo', 0, 'rawHi', 100, ...
                  'fill', 255);                           % % of the value's area from full inversions
VAR.ndsi = struct('dir', fullfile(C.root, 'SnowCover_Daily_TIFFs'), ...
                  'pat', 'NDSI_SnowCover_MOD10A1_daily_%d.tif', 'scale', 0.01, 'rawLo', 0, 'rawHi', 100, ...
                  'fill', -32768);
VAR.ssrd = struct('dir', fullfile(C.root, 'ERA5_Daily_TIFFs'), ...
                  'pat', 'ERA5_ssrd_daily_%d.tif', 'scale', 0.002, 'rawLo', 0, 'rawHi', 25000, ...
                  'fill', -32768);                        % MJ m-2 d-1; 25000 raw = 50 MJ (> any daily sum)

topoTif = fullfile(C.root, 'topography', 'ALOS_AW3D30_topography.tif');
gridOrigin = [-706000 2375000];   % m, upper-left corner of every input raster (EPSG:3338)
gridCell   = 1000;                % m

% ecoregions (settings of alb_ecoregions.m)
E.shp        = fullfile(C.root, 'ecoregions_simplified', 'ecoregions_simplified.shp');
E.field      = 'LEVEL_2';
E.shpCoords  = 'NAD27';        % 'NAD27' (coordinates as in the original coverage) or 'EPSG3338'
E.nad27BBox  = [-2168805.5 410152.1875; 1532912.0 2378649.0];   % [xmin ymin; xmax ymax] of akecoregions
E.datumShift = [-5 135 172];   % m, geocentric NAD27 -> NAD83 shift for Alaska (NIMA TR8350.2)
E.level1 = { ...
    'Alaska Range Transition',      'Boreal'; ...
    'Aleutian Meadows',             'Maritime'; ...
    'Arctic Tundra',                'Polar'; ...
    'Bering Taiga',                 'Polar'; ...
    'Bering Tundra',                'Polar'; ...
    'Coast Mountains Transition',   'Boreal'; ...
    'Coastal Rainforests',          'Maritime'; ...
    'Intermontane Boreal',          'Boreal'; ...
    'Pacific Mountains Transition', 'Boreal'};
E.checkPts = { ...          % (lat, lon WGS84/NAD83) inside a uniform 7 x 7 km block of the region
    'Prudhoe Bay',               70.26, -148.34, 'Arctic Tundra'; ...
    'Fairbanks',                 64.84, -147.72, 'Intermontane Boreal'; ...
    'McGrath',                   62.95, -155.60, 'Intermontane Boreal'; ...
    'Denali summit',             63.07, -151.01, 'Alaska Range Transition'; ...
    'Glennallen',                62.11, -145.55, 'Pacific Mountains Transition'; ...
    'Dillingham',                59.04, -158.46, 'Bering Taiga'; ...
    'Seward Peninsula interior', 65.20, -164.00, 'Bering Tundra'; ...
    'Port Heiden',               56.95, -158.62, 'Aleutian Meadows'; ...
    'Admiralty Island interior', 57.70, -134.30, 'Coastal Rainforests'};
%% =============================================================

assert(C.minPctFull >= 0 && C.minPctFull <= 100, 'C.minPctFull must lie in 0..100');
assert(C.minPctFull == 0 || C.useAlbQA, 'C.minPctFull > 0 needs the QA companion files (C.useAlbQA = true)');
assert(C.minPctFull == 0 || ~strcmpi(C.outRoot, fullfile(C.root, 'Prepared')), ['C.minPctFull > 0 is the ' ...
    'full-inversion sensitivity run: set C.outRoot = fullfile(C.root, ''Prepared_fullinv'') so the primary ' ...
    'layers are not overwritten']);
checkOutRoot(C, VAR);                                  % never overwrite layers made from other albedo files
if ~isfolder(C.outRoot), mkdir(C.outRoot); end
diary(fullfile(C.outRoot, 'prepare_console.txt'));  diary on
cleanDiary = onCleanup(@() diary('off'));
fprintf('alb_prepare started %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
assert(isequal(C.years, unique(C.years)), 'C.years must be sorted and unique');
assert(C.window(1) <= C.opt.pre(1) && C.window(2) >= C.opt.post(2), 'C.window must contain both plateau windows');
assert(C.window(1) >= 60, 'C.window must start on or after 1 Mar (day 60) so that 29 Feb never falls inside');
nY = numel(C.years);  nG = numel(C.aggs);
winDays = (C.window(1):C.window(2))';  nT = numel(winDays);
vUse = {'alb', 'ndsi'};  if C.useERA5, vUse{end+1} = 'ssrd'; end
if C.useAlbQA, vUse{end+1} = 'albqa'; end
for v = vUse                                           % every input file must exist before hours of reading
    miss = C.years(arrayfun(@(y) ~isfile(fullfile(VAR.(v{1}).dir, sprintf(VAR.(v{1}).pat, y))), C.years));
    assert(isempty(miss), 'Missing %s file(s) for %s in %s', v{1}, mat2str(miss), VAR.(v{1}).dir);
end

%% ---------- 1. 1-km reference grid, AOI, topography ----------
refFile = fullfile(VAR.alb.dir, sprintf(VAR.alb.pat, C.years(1)));
meta1 = getMeta(refFile);
assert(~isempty(meta1.geoKeys), 'No GeoKey directory in %s: outputs would carry no CRS', refFile);
assert(isfield(meta1.geoKeys, 'ProjectedCSTypeGeoKey') && meta1.geoKeys.ProjectedCSTypeGeoKey == 3338, ...
    '%s is not tagged EPSG:3338 (ProjectedCSTypeGeoKey)', refFile);
H1 = meta1.H;  W1 = meta1.W;
checkOrigin(refFile, gridOrigin, gridCell);
[topo, Rt] = readgeoraster(topoTif, 'OutputType', 'single');
assert(Rt.CellExtentInWorldX == gridCell && Rt.CellExtentInWorldY == gridCell && ...
       Rt.XWorldLimits(1) == gridOrigin(1) && Rt.YWorldLimits(2) == gridOrigin(2), ...
       '%s is not on the 1-km grid of the albedo files', topoTif);
assert(size(topo, 1) >= H1 && size(topo, 2) >= W1 && size(topo, 3) == 3, '%s: unexpected size', topoTif);
topo = topo(1:H1, 1:W1, :);
aoi1 = isfinite(topo(:,:,1));
assert(~any(any(isfinite(topo(:,:,2)) & ~aoi1)), 'Slope defined outside the AOI');
fprintf('1-km grid %d x %d, AOI %d km2 (slope undefined on %d AOI pixels along coasts and edges)\n', ...
    H1, W1, nnz(aoi1), nnz(aoi1 & ~isfinite(topo(:,:,2))));

%% ---------- 2. ecoregions on the 1-km grid ----------
Sh = shaperead(E.shp);
checkShapeCoords(Sh, E);
recName = cellfun(@strtrim, {Sh.(E.field)}, 'UniformOutput', false);
[regNames, ~, recCode] = unique(recName);    % codes 1..nR, alphabetical
nR = numel(regNames);
regParent = cell(nR, 1);
for r = 1:nR
    k = find(strcmp(E.level1(:,1), regNames{r}));
    assert(numel(k) == 1, 'Region "%s" is not in E.level1', regNames{r});
    regParent{r} = E.level1{k, 2};
end
reg1 = zeros(H1, W1, 'uint8');  nOverlap = 0;
for k = 1:numel(Sh)
    if strcmp(E.shpCoords, 'NAD27')
        [x, y] = nad27AlbersToEPSG3338(Sh(k).X, Sh(k).Y, E.datumShift);
    else
        x = Sh(k).X;  y = Sh(k).Y;
    end
    M = rasterizeEvenOdd(x, y, meta1.R, H1, W1);
    nOverlap = nOverlap + nnz(M & reg1 > 0);
    reg1(M & reg1 == 0) = recCode(k);
end
assert(nOverlap == 0, '%d pixels fall in two ecoregion polygons', nOverlap);
reg1(~aoi1) = 0;
checkLocations(E.checkPts, reg1, regNames, meta1.R);
fprintf('Ecoregions: %d regions; %d km2 of the AOI outside every polygon\n', nR, nnz(aoi1 & reg1 == 0));
clear Sh M x y

%% ---------- 3. aggregated grids and static layers ----------
G = struct([]);
for g = 1:nG
    a = C.aggs(g);  Hc = floor(H1 / a);  Wc = floor(W1 / a);
    dx = a * gridCell;
    Gg.agg = a;  Gg.H = Hc;  Gg.W = Wc;  Gg.minN = ceil(C.minFrac * a^2);
    Gg.dir = fullfile(C.outRoot, sprintf('grid_%dkm', a));
    if ~isfolder(Gg.dir), mkdir(Gg.dir); end
    Gg.meta = meta1;
    Gg.meta.R = maprefcells(gridOrigin(1) + [0, Wc * dx], gridOrigin(2) + [-Hc * dx, 0], [Hc Wc], ...
                            'ColumnsStartFrom', 'north');
    Gg.meta.H = Hc;  Gg.meta.W = Wc;
    blk = @(A) blockSum(A(1:Hc*a, 1:Wc*a), a, Hc, Wc);
    nAoi = blk(aoi1);
    Gg.landFrac = nAoi / a^2;
    Gg.land = nAoi >= ceil(C.landFrac * a^2);
    cnt = zeros(Hc, Wc, nR);
    for r = 1:nR, cnt(:,:,r) = blk(reg1 == r); end
    [mx, arg] = max(cnt, [], 3);
    Gg.reg = zeros(Hc, Wc);  maj = mx >= ceil(C.majorityFrac * a^2);
    Gg.reg(maj) = arg(maj);
    assert(all(Gg.land(Gg.reg > 0)), 'A cell with an ecoregion majority is not a land cell');
    el = topo(:,:,1);  sl = topo(:,:,2);  el(~aoi1) = 0;
    slOk = isfinite(sl);  sl(~slOk) = 0;  nSl = blk(slOk);
    elMean = blk(el) ./ nAoi;
    elSd   = sqrt(max(blk(el.^2) ./ nAoi - elMean.^2, 0) .* nAoi ./ max(nAoi - 1, 1));   % sample SD
    elSd(nAoi < 2) = NaN;
    slMean = blk(sl) ./ nSl;                  % over the AOI pixels with a defined slope
    slMean(nSl == 0) = NaN;
    xc = gridOrigin(1) + dx * ((1:Wc) - 0.5);  yc = gridOrigin(2) - dx * ((1:Hc)' - 0.5);
    [Xc, Yc] = meshgrid(xc, yc);
    [lat, lon] = albersInv(Xc, Yc, [6378137.0, 1/298.257222101], [55 65 50 -154]);   % EPSG:3338 (GRS80)
    Gg.lat = lat;
    Gg.idx = find(Gg.land);                   % land cells, column-major order
    regOut = Gg.reg;  regOut(Gg.reg == 0) = NaN;
    writeTif(fullfile(Gg.dir, 'static_ecoregion.tif'),     regOut, Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_land_frac.tif'),     Gg.landFrac, Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_lat_deg.tif'),       lat, Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_lon_deg.tif'),       lon, Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_elev_mean_m.tif'),   nanOutside(elMean, Gg.land), Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_elev_sd_m.tif'),     nanOutside(elSd, Gg.land), Gg.meta);
    writeTif(fullfile(Gg.dir, 'static_slope_mean_deg.tif'), nanOutside(slMean, Gg.land), Gg.meta);
    codes = struct('code', num2cell(1:nR)', 'level2', regNames(:), 'level1', regParent(:), ...
                   'cells', num2cell(arrayfun(@(r) nnz(Gg.reg == r), 1:nR))', ...
                   'area_km2', num2cell(arrayfun(@(r) nnz(Gg.reg == r) * a^2, 1:nR))');
    writeCsv(fullfile(Gg.dir, 'ecoregion_codes.csv'), codes');
    writeCsv(fullfile(Gg.dir, 'years.csv'), struct('band', num2cell(1:nY), 'year', num2cell(C.years)));
    fprintf('Grid %d km: %d x %d cells, %d land cells, %d with an ecoregion, minN %d\n', ...
        a, Hc, Wc, nnz(Gg.land), nnz(Gg.reg > 0), Gg.minN);
    G = [G, Gg];                                                               %#ok<AGROW>
    clear Gg
end
clear topo el sl cnt mx arg lat lon Xc Yc

%% ---------- 4. per year: read, aggregate, metrics, energy ----------
outNames = {'scd75', 'scd50', 'scd25', 'dur', 'k_rel', 'a_snow', 'a_floor', 'amp', 'nobs', 'qc', ...
            'ndsi_scd50', 'ndsi_qc', 'sw_in_toa', 'sw_abs_toa', 'fill_frac_toa'};
if C.useERA5, outNames = [outNames, {'sw_in_era5', 'sw_abs_era5', 'fill_frac_era5'}]; end
if C.useAlbQA, outNames{end+1} = 'mag_frac'; end
OUT = cell(nG, 1);
for g = 1:nG
    for k = 1:numel(outNames), OUT{g}.(outNames{k}) = nan(numel(G(g).idx), nY, 'single'); end
end
LOG = struct([]);

for iy = 1:nY
    yr = C.years(iy);  t0 = tic;
    qaFile = '';  if C.useAlbQA, qaFile = fullfile(VAR.albqa.dir, sprintf(VAR.albqa.pat, yr)); end
    [SA, NA, doyA, badA, MA, qsA, missA] = aggregateFile(fullfile(VAR.alb.dir, sprintf(VAR.alb.pat, yr)), VAR.alb, C, G, H1, W1, ...
                                                  gridOrigin, gridCell, qaFile, VAR.albqa);
    [SN, NN, doyN, badN, ~, ~, missN] = aggregateFile(fullfile(VAR.ndsi.dir, sprintf(VAR.ndsi.pat, yr)), VAR.ndsi, C, G, H1, W1, ...
                                                       gridOrigin, gridCell);
    assert(isequal(doyA, doyN), '%d: albedo and NDSI windows differ', yr);
    if C.useERA5
        [SS, NS, doyS, badS] = aggregateFile(fullfile(VAR.ssrd.dir, sprintf(VAR.ssrd.pat, yr)), VAR.ssrd, C, G, H1, W1, gridOrigin, gridCell);
        assert(isequal(doyA, doyS), '%d: albedo and SSRD windows differ', yr);
    end
    for g = 1:nG
        id = G(g).idx;
        % ----- albedo -----
        A = SA{g}(id, :) ./ NA{g}(id, :);   A(NA{g}(id, :) < G(g).minN) = NaN;
        v = A(isfinite(A));
        assert(all(v >= 0 & v <= 1), '%d: cell-day albedo outside [0, 1]', yr);
        Mx = transitionMetrics(A', winDays, C.opt);
        qc = reasonCode(Mx, C.opt);
        ok = Mx.ok;
        OUT{g}.scd75(:, iy)   = maskOk(Mx.t75, ok);
        OUT{g}.scd50(:, iy)   = maskOk(Mx.t50, ok);
        OUT{g}.scd25(:, iy)   = maskOk(Mx.t25, ok);
        OUT{g}.dur(:, iy)     = maskOk(Mx.dur, ok);
        OUT{g}.k_rel(:, iy)   = maskOk(0.5 ./ Mx.dur, ok);
        OUT{g}.a_snow(:, iy)  = maskOk(Mx.a_snow, ok);
        OUT{g}.a_floor(:, iy) = maskOk(Mx.a_floor, ok);
        OUT{g}.amp(:, iy)     = maskOk(Mx.amp, ok);
        OUT{g}.nobs(:, iy)    = Mx.nobs(:);
        OUT{g}.qc(:, iy)      = qc(:);
        % ----- NDSI (validation) -----
        B = SN{g}(id, :) ./ NN{g}(id, :);   B(NN{g}(id, :) < G(g).minN) = NaN;
        v = B(isfinite(B));
        assert(all(v >= 0 & v <= 1), '%d: cell-day NDSI snow cover outside [0, 1]', yr);
        Mn = transitionMetrics(B', winDays, C.opt);
        OUT{g}.ndsi_scd50(:, iy) = maskOk(Mn.t50, Mn.ok);
        OUT{g}.ndsi_qc(:, iy)    = reshape(reasonCode(Mn, C.opt), [], 1);
        % ----- shortwave -----
        Q = toaDaily(G(g).lat(id)', doyA(:), C.S0) / 1e6;             % nT x nLand, MJ m-2 d-1
        [swIn, swAbs, ff] = springEnergy(Mx, Q);
        OUT{g}.sw_in_toa(:, iy) = swIn(:);  OUT{g}.sw_abs_toa(:, iy) = swAbs(:);  OUT{g}.fill_frac_toa(:, iy) = ff(:);
        nE5 = NaN;
        if C.useERA5
            SW = SS{g}(id, :) ./ NS{g}(id, :);   SW(NS{g}(id, :) < 1) = NaN;
            SW = SW';
            v = SW(isfinite(SW));
            assert(all(v >= 0 & v <= 50), '%d: daily SSRD outside [0, 50] MJ m-2', yr);
            [swIn, swAbs, ff] = springEnergy(Mx, SW);
            OUT{g}.sw_in_era5(:, iy) = swIn(:);  OUT{g}.sw_abs_era5(:, iy) = swAbs(:);  OUT{g}.fill_frac_era5(:, iy) = ff(:);
            nE5 = nnz(isfinite(swAbs));
            assert(all(swIn(isfinite(swIn)) <= sum(Q(:, isfinite(swIn)), 1) + 1e-6), ...
                '%d: ERA5 surface shortwave exceeds the TOA insolation', yr);
        end
        % ----- magnitude-inversion share of the observed area -----
        magPct = NaN;
        if C.useAlbQA
            [mf, magPct] = magFraction(NA{g}(id, :), MA{g}(id, :), G(g).minN);
            v = mf(isfinite(mf));
            assert(all(v >= 0 & v <= 1 + 1e-12), '%d: mag_frac outside [0, 1]', yr);
            OUT{g}.mag_frac(:, iy) = mf;
        end
        row = struct('year', yr, 'grid_km', G(g).agg, 'land_cells', numel(id), ...
            'alb_celldays_valid_pct', 100 * mean(isfinite(A(:))), 'alb_ok_cells', nnz(ok), ...
            'alb_ok_pct', 100 * nnz(ok) / numel(id), 'qc1_few_days', nnz(qc == 1), 'qc2_no_plateau', nnz(qc == 2), ...
            'qc3_no_snow_signal', nnz(qc == 3), 'qc4_no_crossing', nnz(qc == 4), ...
            'ndsi_ok_cells', nnz(Mn.ok), 'energy_toa_cells', nnz(isfinite(OUT{g}.sw_abs_toa(:, iy))), ...
            'energy_era5_cells', nE5, 'alb_raw_out_of_range', badA, 'ndsi_raw_out_of_range', badN, ...
            'alb_mag_area_pct', magPct, 'alb_qa_dropped', qsA.dropped, 'alb_qa_mismatch', qsA.mismatch, ...
            'alb_window_days_missing', numel(missA), 'ndsi_window_days_missing', numel(missN), ...
            'seconds', toc(t0));
        LOG = [LOG, row];                                                    %#ok<AGROW>
        fprintf('  %d  %2d km: albedo transition valid in %5.1f %% of land cells (NDSI %5.1f %%)\n', ...
            yr, G(g).agg, row.alb_ok_pct, 100 * nnz(Mn.ok) / numel(id));
        if C.useAlbQA
            fprintf('            magnitude inversions: %4.1f %% of the observed area in kept cell-days\n', magPct);
        end
    end
    if badA > 0 || badN > 0
        fprintf('  %d: raw values outside the valid range and not fill: albedo %d, NDSI %d (treated as missing)\n', yr, badA, badN);
    end
    if C.useAlbQA && (qsA.dropped > 0 || qsA.mismatch > 0)
        fprintf('  %d: albedo values dropped by C.minPctFull %d; valid albedo values without a QA value %d\n', ...
            yr, qsA.dropped, qsA.mismatch);
    end
    fprintf('%d done (%.0f s)\n', yr, toc(t0));
    clear SA NA MA SN NN SS NS A B SW Mx Mn Q
end

%% ---------- 5. write ----------
for g = 1:nG
    Hc = G(g).H;  Wc = G(g).W;
    for k = 1:numel(outNames)
        A = nan(Hc * Wc, nY, 'single');
        A(G(g).idx, :) = OUT{g}.(outNames{k});
        writeTif(fullfile(G(g).dir, [outNames{k} '.tif']), reshape(A, Hc, Wc, nY), G(g).meta);
    end
    writeCsv(fullfile(G(g).dir, 'prepare_log.csv'), LOG([LOG.grid_km] == G(g).agg));
    writeSettings(fullfile(G(g).dir, 'settings.csv'), C, G(g), VAR, outNames);
    fprintf('Wrote %d layers to %s\n', numel(outNames), G(g).dir);
end
fprintf('alb_prepare finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));


%% =====================================================================
%%                     READING AND AGGREGATION
%% =====================================================================

function [S, N, doyAct, nBad, M, qs, missing] = aggregateFile(fp, V, C, G, H1, W1, origin, cellSize, qaFp, QV)
% One yearly GeoTIFF -> block sums S{g} and valid-pixel counts N{g} per grid g,
% [cells x window days] (cells in column-major order of grid g). Read one tile
% (all bands) at a time; only the window bands are kept.
% Albedo with its QA companion file qaFp (QV = VAR.albqa): values whose
% full-inversion share is below C.minPctFull are dropped (qaFilter), and M{g}
% sums the magnitude-inversion share (1 - pct_full/100) of the kept values.
% qs.dropped / qs.mismatch / qs.valid count the values dropped by that rule,
% the valid values without a QA value, and the kept values (NaN without qaFp).
if nargin < 9, qaFp = ''; end
hasQa = ~isempty(qaFp);
assert(isfile(fp), 'Missing input file %s', fp);
info = imfinfo(fp);  info = info(1);
checkOrigin(fp, origin, cellSize);
assert(info.Height >= H1 && info.Width >= W1, '%s is smaller than the albedo grid', fp);
nB = info.SamplesPerPixel;
yr = str2double(regexp(fp, '(\d{4})\.tif$', 'tokens', 'once'));
dn = bandDates(info, fp, nB, yr, C.bandDates);
[bIdx, doyAct, pos, missing] = windowBands(dn, C.window, fp);
nGap = round(dn(end) - dn(1)) + 1 - numel(dn);     % days without a band between the first and last band
if nGap > 0 || ~isempty(missing)
    ms = 'none';  if ~isempty(missing), ms = mat2str(missing); end
    fprintf('  note: %s lacks %d day(s) between %s and %s; window days without a band: %s (no observations)\n', ...
        fp, nGap, datestr(dn(1), 'yyyy-mm-dd'), datestr(dn(end), 'yyyy-mm-dd'), ms);
end
nT = numel(doyAct);                                 % full window; days without a band stay 0 / 0
t = Tiff(fp, 'r');  closer = onCleanup(@() close(t));
Fk = struct('t', t, 'fp', fp, 'nB', nB);
tw = getTag(t, 'TileWidth');  th = getTag(t, 'TileLength');
qs = struct('dropped', NaN, 'mismatch', NaN, 'valid', NaN);
if hasQa
    assert(isfile(qaFp), 'Missing QA file %s', qaFp);
    qinfo = imfinfo(qaFp);  qinfo = qinfo(1);
    checkOrigin(qaFp, origin, cellSize);
    assert(qinfo.Height == info.Height && qinfo.Width == info.Width, '%s and %s differ in size', qaFp, fp);
    qn = bandDates(qinfo, qaFp, qinfo.SamplesPerPixel, yr, C.bandDates);
    % The QA file must hold every window band the albedo file uses; QA columns are
    % matched to those bands by date. Dates outside them may differ (a partial year
    % whose two files were exported on different days: GEE keeps adding MCD43A3 days).
    [inQ, qIdx] = ismember(dn(bIdx), qn);
    assert(all(inQ), '%s lacks window band(s) %s that %s has', qaFp, ...
        strjoin(cellstr(datestr(dn(bIdx(~inQ)), 'yyyy-mm-dd')), ', '), fp);
    if ~isequal(qn, dn)
        dd = setxor(qn, dn);
        fprintf('  note: %s and %s differ in %d band date(s) (%s ... %s), none of them a window band used: ignored\n', ...
            qaFp, fp, numel(dd), datestr(dd(1), 'yyyy-mm-dd'), datestr(dd(end), 'yyyy-mm-dd'));
    end
    tq = Tiff(qaFp, 'r');  closerQ = onCleanup(@() close(tq));
    Fq = struct('t', tq, 'fp', qaFp, 'nB', qinfo.SamplesPerPixel);
    try
        sameTiles = getTag(tq, 'TileWidth') == tw && getTag(tq, 'TileLength') == th;
    catch
        sameTiles = false;                              % QA file not tiled
    end
    if ~sameTiles
        fprintf('  note: %s is tiled differently from the albedo file; reading it by pixel region (slower)\n', qaFp);
    end
    qs = struct('dropped', 0, 'mismatch', 0, 'valid', 0);
end
nG = numel(G);  S = cell(nG, 1);  N = cell(nG, 1);  M = cell(nG, 1);
for g = 1:nG
    S{g} = zeros(G(g).H * G(g).W, nT);  N{g} = zeros(G(g).H * G(g).W, nT);
    if hasQa, M{g} = zeros(G(g).H * G(g).W, nT); end
end
Hmax = max([G.H] .* [G.agg]);  Wmax = max([G.W] .* [G.agg]);   % rows / columns used by any grid
nBad = 0;
for ty = 0:ceil(Hmax / th) - 1
    rr = ty*th + 1 : min((ty+1)*th, info.Height);
    for tx = 0:ceil(Wmax / tw) - 1
        cc = tx*tw + 1 : min((tx+1)*tw, info.Width);
        a = readTile(Fk, rr, cc);                  % pixels x bands, int16
        a = a(:, bIdx);
        m = a >= V.rawLo & a <= V.rawHi;
        isFill = false(size(a));
        for fv = V.fill, isFill = isFill | a == fv; end
        nBad = nBad + nnz(~m & ~isFill);
        if hasQa
            if sameTiles
                q = readTile(Fq, rr, cc);
            else
                q = imread(qaFp, 'PixelRegion', {[rr(1) rr(end)], [cc(1) cc(end)]});
                assert(size(q, 3) == Fq.nB, 'Band count mismatch reading %s', qaFp);
                q = reshape(q, [], size(q, 3));
            end
            [m, mag, nDrop, nMis] = qaFilter(m, double(q(:, qIdx)), QV, C.minPctFull);   % QA bands matched by date
            qs.dropped = qs.dropped + nDrop;  qs.mismatch = qs.mismatch + nMis;  qs.valid = qs.valid + nnz(m);
        end
        x = double(a) * V.scale;  x(~m) = 0;
        [RI, CI] = ndgrid(rr, cc);  RI = RI(:);  CI = CI(:);
        for g = 1:nG
            ag = G(g).agg;
            in = RI <= G(g).H * ag & CI <= G(g).W * ag;
            if ~any(in), continue; end
            lin = floor((CI(in) - 1) / ag) * G(g).H + floor((RI(in) - 1) / ag) + 1;
            [u, ~, j] = unique(lin);
            P = sparse(j, 1:numel(j), 1, numel(u), numel(j));
            S{g}(u, pos) = S{g}(u, pos) + P * x(in, :);
            N{g}(u, pos) = N{g}(u, pos) + P * double(m(in, :));
            if hasQa, M{g}(u, pos) = M{g}(u, pos) + P * mag(in, :); end
        end
    end
end
for g = 1:nG
    assert(max(N{g}(:)) <= G(g).agg^2, '%s: more valid pixels than a cell holds', fp);
end
if hasQa       % the two files come from the same GEE mask: disagreement means a wrong pairing
    assert(qs.mismatch <= 1e-3 * max(qs.valid, 1), ['%s: %d valid albedo values have no QA value ' ...
        '(%.2f %% of the kept values) - is %s the companion of this file?'], fp, qs.mismatch, ...
        100 * qs.mismatch / max(qs.valid, 1), qaFp);
end
end

function [m, mag, nDrop, nMis] = qaFilter(m, q, QV, minPct)
% m: valid albedo values (P x T logical); q: QA companion values, % of each
% value's area that comes from full BRDF inversions (QV.rawLo..QV.rawHi valid,
% anything else = no QA value). Keeps the values with q >= minPct (minPct = 0:
% all values, including any without a QA value). mag = magnitude-inversion
% share 1 - q/100 of each kept value with a QA value (0 elsewhere).
% nDrop = valid values with a QA value below minPct; nMis = valid values
% without a QA value (also dropped when minPct > 0).
qOk = q >= QV.rawLo & q <= QV.rawHi;
nMis = nnz(m & ~qOk);
nDrop = 0;
if minPct > 0
    low = m & qOk & q < minPct;
    nDrop = nnz(low);
    m = m & qOk & ~low;
end
mag = zeros(size(q));
k = m & qOk;
mag(k) = 1 - q(k) / 100;
end

function [mf, pct] = magFraction(N, M, minN)
% Magnitude-inversion share of the observed area per cell (rows), over the
% cell-days kept for the albedo series (N >= minN): sum(M) / sum(N) on those
% days; NaN where no cell-day is kept. pct = the same share over all cells, %.
use = N >= minN;
n = sum(N .* use, 2);  s = sum(M .* use, 2);
mf = s ./ n;  mf(n == 0) = NaN;
pct = NaN;
if sum(n) > 0, pct = 100 * sum(s) / sum(n); end
end

function a = readTile(Fk, rr, cc)
% One tile, all bands -> P x nBands int16. Falls back to imread if needed.
try
    a = readEncodedTile(Fk.t, computeTile(Fk.t, [rr(1) cc(1)]));
    a = a(1:numel(rr), 1:numel(cc), :);
catch
    a = imread(Fk.fp, 'PixelRegion', {[rr(1) rr(end)], [cc(1) cc(end)]});
end
assert(size(a,3) == Fk.nB, 'Band count mismatch reading %s', Fk.fp);
a = reshape(a, [], size(a,3));
end

function dn = bandDates(info, fp, nB, yr, mode)
% Band dates from the GDAL band descriptions 'dYYYYMMDD' (TIFF tag 42112), or
% (mode 'dayOfYear') band k = day k of year yr. Stops if the dates are missing,
% incomplete, not strictly increasing or not in year yr. Days without a band
% are allowed; windowBands treats those inside the window as days without
% observations.
dn = [];
if strcmp(mode, 'dayOfYear')
    dn = datenum(yr, 1, 1) + (0:nB-1)';
    assert(nB == datenum(yr + 1, 1, 1) - datenum(yr, 1, 1), '%s: %d bands but %d has %d days', ...
        fp, nB, yr, datenum(yr + 1, 1, 1) - datenum(yr, 1, 1));
    return
end
assert(strcmp(mode, 'tag'), 'C.bandDates must be ''tag'' or ''dayOfYear''');
% The GDAL_METADATA tag is read from the file itself (gdalMetadata): imfinfo
% does not return it for BigTIFF files, and every GEE export here is BigTIFF.
tok = regexp(gdalMetadata(fp), 'sample="(\d+)"[^>]*role="description">d(\d{8})<', 'tokens');
if ~isempty(tok)
    smp = cellfun(@(c) str2double(c{1}), tok);  v = cellfun(@(c) str2double(c{2}), tok);
    [smp, k] = sort(smp);  v = v(k);
    assert(isequal(smp, 0:numel(smp) - 1), '%s: band descriptions do not cover bands 1..%d', fp, numel(smp));
    dn = datenum(floor(v/1e4), mod(floor(v/100), 100), mod(v, 100));
end
dn = dn(:);
assert(numel(dn) == nB, ['%s: %d band dates found for %d bands (GDAL band descriptions dYYYYMMDD ' ...
    'expected). If the files have one band per day from 1 Jan, set C.bandDates = ''dayOfYear''.'], fp, numel(dn), nB);
assert(all(diff(dn) > 0), '%s: band dates are not strictly increasing', fp);   % gaps: see windowBands
dv = datevec(dn);
assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fp, yr);
end

function [bIdx, doyAct, pos, missing] = windowBands(dn, win, fp)
% Bands whose date falls on days win(1)..win(2) of a 365-day calendar (leap
% years: 29 Feb skipped, later days shifted by one), and their columns pos in
% the full window. Window days without a band (MODIS outages) are listed in
% missing and enter the series as days without observations. doyAct = actual
% day of year of every window day (solar geometry), from the calendar, so all
% files of a year share it.
assert(win(1) >= 60, 'windowBands: the window must start on or after day 60 (1 Mar)');
dv = datevec(dn);
yr = dv(1, 1);
leap = mod(yr, 4) == 0 && (mod(yr, 100) ~= 0 || mod(yr, 400) == 0);
doy365 = dn - datenum(yr, 1, 0) - (leap & dv(:, 2) >= 3);
isFeb29 = dv(:, 2) == 2 & dv(:, 3) == 29;
bIdx = find(doy365 >= win(1) & doy365 <= win(2) & ~isFeb29);
assert(~isempty(bIdx), '%s: no band falls inside the window', fp);
wd = (win(1):win(2))';
[~, pos] = ismember(doy365(bIdx), wd);
missing = setdiff(wd, doy365(bIdx))';
doyAct = wd + leap;                                % every window day is after 29 Feb
end

function s = gdalMetadata(fp)
% Text of the GDAL_METADATA tag (42112) of the first image directory, read
% straight from the file because imfinfo does not return it for BigTIFF files
% (the GEE export format). Classic TIFF and BigTIFF, either byte order.
% '' if the tag is absent.
fid = fopen(fp, 'r');  assert(fid > 0, 'Cannot open %s', fp);
closer = onCleanup(@() fclose(fid));
bo = fread(fid, 2, 'uint8=>char')';
if strcmp(bo, 'II'), mf = 'ieee-le'; elseif strcmp(bo, 'MM'), mf = 'ieee-be';
else, error('%s is not a TIFF file', fp); end
ver = fread(fid, 1, 'uint16', 0, mf);
if ver == 42                                         % classic TIFF: 4-byte counts and offsets
    cntT = 'uint32';  offT = 'uint32';  nT = 'uint16';  entry = 12;  inl = 4;
    ifd = fread(fid, 1, offT, 0, mf);
elseif ver == 43                                     % BigTIFF: 8-byte counts and offsets
    cntT = 'uint64';  offT = 'uint64';  nT = 'uint64';  entry = 20;  inl = 8;
    fread(fid, 2, 'uint16', 0, mf);                  % offset size (8), reserved (0)
    ifd = fread(fid, 1, offT, 0, mf);
else
    error('%s: unknown TIFF version %d', fp, ver);
end
fseek(fid, ifd, 'bof');
n = fread(fid, 1, nT, 0, mf);
s = '';
for k = 1:n
    e = ftell(fid);
    tag = fread(fid, 1, 'uint16', 0, mf);  typ = fread(fid, 1, 'uint16', 0, mf);
    cnt = fread(fid, 1, cntT, 0, mf);
    if tag == 42112 && (typ == 2 || typ == 7)
        if cnt > inl, fseek(fid, fread(fid, 1, offT, 0, mf), 'bof'); end
        s = fread(fid, cnt, 'uint8=>char')';
        break
    end
    fseek(fid, e + entry, 'bof');
end
s = s(s ~= char(0));
end

function checkOrigin(fp, origin, cellSize)
% Upper-left corner and cell size of the GeoTIFF must match the study grid.
gi = georasterinfo(fp);  R = gi.RasterReference;
x0 = R.XWorldLimits(1);  y0 = R.YWorldLimits(2);
assert(x0 == origin(1) && y0 == origin(2) && R.CellExtentInWorldX == cellSize && R.CellExtentInWorldY == cellSize, ...
    '%s: grid origin (%g, %g) / cell %g m differ from (%g, %g) / %g m', fp, x0, y0, R.CellExtentInWorldX, ...
    origin(1), origin(2), cellSize);
end

function A = nanOutside(A, keep)
A(~keep) = NaN;
end

function S = blockSum(B, a, H, W)
% Sum of B over non-overlapping a x a blocks starting at the top-left corner.
S = reshape(sum(sum(reshape(double(B), a, H, a, W), 1), 3), H, W);
end


%% =====================================================================
%%                  SPRING TRANSITION AND SHORTWAVE
%% =====================================================================

function M = transitionMetrics(A, doy, opt)
% Spring albedo-transition metrics (numerics of alb_metrics.m; every option
% must be given). A: [nT x nC] one year, NaN = no retrieval; doy: [nT x 1].
%   Step 1  Gaussian kernel regression on the irregularly observed days
%           (Nadaraya-Watson), as two matrix products.
%   Step 2  Plateaux from robust percentiles of the smoothed curve, then
%           threshold crossings.
need = {'bw', 'minwt', 'pre', 'post', 'qhi', 'qlo', 'minamp', 'minobs'};
assert(all(isfield(opt, need)), 'transitionMetrics: options missing');
doy = doy(:);  nT = numel(doy);  nC = size(A, 2);
assert(size(A, 1) == nT, 'A must be nT x nC with nT = numel(doy)');

% ---- 1. kernel smoothing ----
K  = exp(-0.5*((doy - doy.')/opt.bw).^2);        % nT x nT
Vm = ~isnan(A);
Az = A;  Az(~Vm) = 0;
num = K * Az;
den = K * double(Vm);
S   = num ./ den;
S(den < opt.minwt) = NaN;

% ---- 2. plateaux and crossings ----
ipre  = doy >= opt.pre(1)  & doy <= opt.pre(2);
ipost = doy >= opt.post(1) & doy <= opt.post(2);
a_snow  = nanpctl(S(ipre ,:), opt.qhi);
a_floor = nanpctl(S(ipost,:), opt.qlo);
amp     = a_snow - a_floor;

t75 = cross_level(S, doy, a_floor + 0.75*amp);
t50 = cross_level(S, doy, a_floor + 0.50*amp);
t25 = cross_level(S, doy, a_floor + 0.25*amp);
dur = t25 - t75;

nobs = sum(Vm, 1);
ok = isfinite(t50) & isfinite(dur) & amp > opt.minamp & dur > 0 & dur < 80 ...
     & nobs >= opt.minobs & isfinite(a_snow) & isfinite(a_floor);

M = struct();
M.a_snow = a_snow;  M.a_floor = a_floor;  M.amp = amp;
M.t50 = t50;        M.t75 = t75;          M.t25 = t25;
M.dur = dur;        M.nobs = nobs;
M.ok  = ok;         M.S = S;              M.doy = doy;
end

function t = cross_level(S, doy, lev)
% First downward crossing of LEV, linearly interpolated between the bracketing
% fitted days.  lev is 1 x nC.
[nT,nC] = size(S);
below = S < lev & ~isnan(S);
[has,k] = max(below,[],1);              % k = first TRUE row
t = nan(1,nC);
j = find(has & k>1);
if isempty(j), return; end
lin1 = sub2ind([nT nC], k(j)-1, j);
lin2 = sub2ind([nT nC], k(j)  , j);
y1 = S(lin1);  y2 = S(lin2);
t1 = doy(k(j)-1).';  t2 = doy(k(j)).';
w = (y1 - lev(j)) ./ (y1 - y2);
w(~isfinite(w)) = 0.5;
t(j) = t1 + min(max(w,0),1) .* (t2 - t1);
end

function q = nanpctl(X, p)
% Column-wise percentile ignoring NaN.  No Statistics Toolbox required.
[nR,nC] = size(X);
q = nan(1,nC);
for c = 1:nC
    v = sort(X(~isnan(X(:,c)),c));
    n = numel(v);
    if n==0, continue; end
    if n==1, q(c) = v; continue; end
    % MATLAB's prctile convention: midpoints of the empirical cdf
    pos = (p/100)*n - 0.5 + 1;
    pos = min(max(pos,1),n);
    lo  = floor(pos);  hi = ceil(pos);
    q(c) = v(lo) + (pos-lo)*(v(hi)-v(lo));
end
end

function w = reasonCode(M, opt)
% QC outcome per cell: 0 valid; 1 fewer than opt.minobs valid days; 2 a plateau
% window has no smoothed value; 3 amplitude <= opt.minamp (no seasonal snow
% signal); 4 no crossing, or duration outside (0, 80) d. Earlier causes win.
w = zeros(size(M.ok), 'uint8');
bad = ~M.ok;
w(bad) = 4;
w(bad & ~(M.amp > opt.minamp)) = 3;
w(bad & ~(isfinite(M.a_snow) & isfinite(M.a_floor))) = 2;
w(bad & M.nobs < opt.minobs) = 1;
end

function v = maskOk(x, ok)
% Column vector of x with NaN where the cell-year failed QC.
v = x(:);  v(~ok(:)) = NaN;
end

function [swIn, swAbs, fillFrac] = springEnergy(M, SW)
% Window totals of incoming and absorbed shortwave for QC-passed cell-years.
% M: transitionMetrics output (S = smoothed albedo, nT x nC); SW: nT x nC daily
% incoming shortwave, MJ m-2 d-1. A cell-year also needs SW on every day.
% Gaps in S are filled: inside the record linearly in time; before the first
% smoothed day with a_snow; after the last with a_floor. fillFrac = share of
% sum(SW) on filled days.
[nT, nC] = size(M.S);
assert(isequal(size(SW), [nT nC]), 'springEnergy: SW must be nT x nC');
swIn = nan(1, nC);  swAbs = nan(1, nC);  fillFrac = nan(1, nC);
j = find(M.ok & all(isfinite(SW), 1));
if isempty(j), return; end
S = M.S(:, j);  F = isnan(S);  t = (1:nT)';
for k = find(any(F, 1))
    s = S(:, k);  v = ~F(:, k);
    i1 = find(v, 1, 'first');  i2 = find(v, 1, 'last');
    s(1:i1-1) = M.a_snow(j(k));  s(i2+1:end) = M.a_floor(j(k));
    gap = F(:, k);  gap(1:i1-1) = false;  gap(i2+1:end) = false;
    if any(gap), s(gap) = interp1(t(v), s(v), t(gap), 'linear'); end
    S(:, k) = s;
end
assert(all(isfinite(S(:))) && all(S(:) >= 0 & S(:) <= 1), 'springEnergy: filled albedo outside [0, 1]');
w = SW(:, j);
assert(all(w(:) >= 0), 'springEnergy: negative shortwave');
swIn(j)  = sum(w, 1);
assert(all(swIn(j) > 0), 'springEnergy: no incoming shortwave in the window');
swAbs(j) = sum(w .* (1 - S), 1);
fillFrac(j) = sum(w .* F, 1) ./ swIn(j);
end

function Q = toaDaily(latDeg, doy, S0)
% Daily top-of-atmosphere insolation, J m-2 d-1: Spencer (1971) declination
% and eccentricity factor (Iqbal 1983, eqs 1.2.1 and 1.3.1), sunset hour angle
% clipped for polar day and polar night. doy may be fractional (non-leap year).
G   = 2 * pi * (doy - 1) / 365;
dec = 0.006918 - 0.399912 * cos(G) + 0.070257 * sin(G) - 0.006758 * cos(2*G) + 0.000907 * sin(2*G) ...
      - 0.002697 * cos(3*G) + 0.00148 * sin(3*G);
E0  = 1.000110 + 0.034221 * cos(G) + 0.001280 * sin(G) + 0.000719 * cos(2*G) + 0.000077 * sin(2*G);
phi = latDeg * pi / 180;
ws  = acos(min(max(-tan(phi) .* tan(dec), -1), 1));
Q   = 86400 / pi * S0 * E0 .* (ws .* sin(phi) .* sin(dec) + cos(phi) .* cos(dec) .* sin(ws));
end


%% =====================================================================
%%            ECOREGION GEOMETRY (verbatim from alb_ecoregions.m)
%% =====================================================================

function checkShapeCoords(Sh, E)
x = [Sh.X];  y = [Sh.Y];
bb = [min(x) min(y); max(x) max(y)];
switch E.shpCoords
    case 'NAD27'
        assert(max(abs(bb(:) - E.nad27BBox(:))) < 1, ['The shapefile extent differs from the original NAD27 ' ...
            'akecoregions extent by %.1f m: the coordinates may have been reprojected. Check the file and ' ...
            'set E.shpCoords = ''EPSG3338'' if they are true NAD83 Alaska Albers.'], max(abs(bb(:) - E.nad27BBox(:))));
    case 'EPSG3338'
    otherwise
        error('E.shpCoords must be ''NAD27'' or ''EPSG3338''');
end
end

function [x2, y2] = nad27AlbersToEPSG3338(x, y, shift)
% NAD27 / Alaska Albers (Clarke 1866) metres -> NAD83 / Alaska Albers (EPSG:3338).
% Albers equal-area formulas: Snyder (1987), USGS PP 1395, eqs. 14-1..14-21.
clk = [6378206.4, 1 - 6356583.8/6378206.4];    % Clarke 1866: a, f
grs = [6378137.0, 1/298.257222101];            % GRS80
alb = [55 65 50 -154];                         % phi1 phi2 phi0 lambda0 (deg), FE = FN = 0
nn = ~isnan(x);
x2 = nan(size(x));  y2 = nan(size(y));
[lat, lon] = albersInv(x(nn), y(nn), clk, alb);
[X, Y, Z]  = geod2ecef(lat, lon, clk);
[lat, lon] = ecef2geod(X + shift(1), Y + shift(2), Z + shift(3), grs);
[x2(nn), y2(nn)] = albersFwd(lat, lon, grs, alb);
end

function [n, C, r0, e] = albersConst(ell, alb)
a = ell(1);  f = ell(2);  e = sqrt(2*f - f^2);  d = pi/180;
q = @(p) (1 - e^2) * (sin(p) ./ (1 - e^2*sin(p).^2) - 1/(2*e) * log((1 - e*sin(p)) ./ (1 + e*sin(p))));
mm = @(p) cos(p) ./ sqrt(1 - e^2*sin(p).^2);
p1 = alb(1)*d;  p2 = alb(2)*d;  p0 = alb(3)*d;
n  = (mm(p1)^2 - mm(p2)^2) / (q(p2) - q(p1));
C  = mm(p1)^2 + n*q(p1);
r0 = a * sqrt(C - n*q(p0)) / n;
end

function [x, y] = albersFwd(lat, lon, ell, alb)
[n, C, r0, e] = albersConst(ell, alb);  a = ell(1);  d = pi/180;
s = sin(lat*d);
q = (1 - e^2) * (s ./ (1 - e^2*s.^2) - 1/(2*e) * log((1 - e*s) ./ (1 + e*s)));
r = a * sqrt(C - n*q) / n;  th = n * (lon - alb(4)) * d;
x = r .* sin(th);  y = r0 - r .* cos(th);
end

function [lat, lon] = albersInv(x, y, ell, alb)
[n, C, r0, e] = albersConst(ell, alb);  a = ell(1);  d = pi/180;
r = hypot(x, r0 - y);  th = atan2(x, r0 - y);
q = (C - (r*n/a).^2) / n;
phi = asin(q/2);
for it = 1:15                                   % Snyder eq. 14-19 (converges in <5 iterations)
    s = sin(phi);
    dphi = (1 - e^2*s.^2).^2 ./ (2*cos(phi)) .* (q/(1 - e^2) - s./(1 - e^2*s.^2) + 1/(2*e)*log((1 - e*s)./(1 + e*s)));
    phi = phi + dphi;
end
assert(max(abs(dphi(:))) < 1e-12, 'Albers inverse did not converge');   % dphi is a matrix for grid input
lat = phi / d;  lon = alb(4) + th / (n*d);
end

function [X, Y, Z] = geod2ecef(lat, lon, ell)
a = ell(1);  f = ell(2);  e2 = 2*f - f^2;  d = pi/180;
N = a ./ sqrt(1 - e2*sin(lat*d).^2);
X = N .* cos(lat*d) .* cos(lon*d);  Y = N .* cos(lat*d) .* sin(lon*d);  Z = N*(1 - e2) .* sin(lat*d);
end

function [lat, lon] = ecef2geod(X, Y, Z, ell)
a = ell(1);  f = ell(2);  e2 = 2*f - f^2;
lon = atan2(Y, X);  p = hypot(X, Y);  lat = atan2(Z, p*(1 - e2));
for it = 1:10
    N = a ./ sqrt(1 - e2*sin(lat).^2);  h = p ./ cos(lat) - N;
    lat = atan2(Z, p .* (1 - e2*N./(N + h)));
end
lat = lat * 180/pi;  lon = lon * 180/pi;
end

function M = rasterizeEvenOdd(x, y, R, H, W)
% Polygon (NaN-separated rings, EPSG:3338 m) -> logical H x W mask.
% Pixel-centre rule; even-odd across all rings (holes are removed).
u = (x(:) - R.XWorldLimits(1)) / R.CellExtentInWorldX;    % 0 = left edge of column 1
v = (R.YWorldLimits(2) - y(:)) / R.CellExtentInWorldY;    % 0 = top edge of row 1
cut = [0; find(isnan(u)); numel(u) + 1];
U1 = []; U2 = []; V1 = []; V2 = [];
for k = 1:numel(cut) - 1
    i = cut(k) + 1 : cut(k+1) - 1;
    if numel(i) < 3, continue; end
    ur = u(i);  vr = v(i);
    ur(end+1) = ur(1);  vr(end+1) = vr(1);                %#ok<AGROW>  close the ring
    U1 = [U1; ur(1:end-1)];  U2 = [U2; ur(2:end)];        %#ok<AGROW>
    V1 = [V1; vr(1:end-1)];  V2 = [V2; vr(2:end)];        %#ok<AGROW>
end
M = false(H, W);
if isempty(V1), return; end
for r = max(1, floor(min(V1)) + 1) : min(H, ceil(max(V1)))
    vc = r - 0.5;                                          % centre line of row r
    k  = (V1 <= vc) ~= (V2 <= vc);                         % edges crossing it (half-open)
    if ~any(k), continue; end
    xi = sort(U1(k) + (vc - V1(k)) .* (U2(k) - U1(k)) ./ (V2(k) - V1(k)));
    assert(mod(numel(xi), 2) == 0, 'Odd number of ring crossings on row %d', r);
    c0 = max(1, ceil(xi(1:2:end) + 0.5));                  % first centre >= left crossing
    c1 = min(W, ceil(xi(2:2:end) + 0.5) - 1);              % last centre <  right crossing
    for j = 1:numel(c0)
        if c1(j) >= c0(j), M(r, c0(j):c1(j)) = true; end
    end
end
end

function checkLocations(pts, reg, regNames, R)
for k = 1:size(pts, 1)
    [x, y] = albersFwd(pts{k,2}, pts{k,3}, [6378137.0, 1/298.257222101], [55 65 50 -154]);
    c = floor((x - R.XWorldLimits(1)) / R.CellExtentInWorldX) + 1;
    r = floor((R.YWorldLimits(2) - y) / R.CellExtentInWorldY) + 1;
    code = reg(r, c);
    got = '(none)';  if code > 0, got = regNames{code}; end
    assert(strcmp(got, pts{k,4}), 'Location check: %s falls in %s, expected %s', pts{k,1}, got, pts{k,4});
end
fprintf('Location checks passed (%d places)\n', size(pts, 1));
end


%% =====================================================================
%%                              OUTPUT
%% =====================================================================

function meta = getMeta(fp)
gi     = georasterinfo(fp);
meta.R = gi.RasterReference;
meta.H = gi.RasterSize(1);
meta.W = gi.RasterSize(2);
try
    gt = geotiffinfo(fp);
    meta.geoKeys = gt.GeoTIFFTags.GeoKeyDirectoryTag;
catch
    meta.geoKeys = [];
end
% geotiffinfo names GeoKeys it does not recognise 'Unknown', and geotiffwrite
% refuses them. GEE exports carry one ([0 0 0], most likely the NAD83 -> WGS84
% TOWGS84 shift); ProjectedCSTypeGeoKey = 3338 defines the CRS on its own.
if isstruct(meta.geoKeys) && isfield(meta.geoKeys, 'Unknown')
    meta.geoKeys = rmfield(meta.geoKeys, 'Unknown');
end
end

function writeTif(fp, A, meta)
A = single(A);
if ~isempty(meta.geoKeys)
    geotiffwrite(fp, A, meta.R, 'GeoKeyDirectoryTag', meta.geoKeys);
else
    geotiffwrite(fp, A, meta.R);
end
end

function writeCsv(fp, S)
% Struct array -> CSV (numeric fields %.8g, char fields quoted).
f = fieldnames(S);
fid = fopen(fp, 'w');  assert(fid > 0, 'Cannot write %s', fp);
fprintf(fid, '%s\n', strjoin(f', ','));
for k = 1:numel(S)
    c = cell(1, numel(f));
    for j = 1:numel(f)
        v = S(k).(f{j});
        if ischar(v), c{j} = ['"' v '"'];
        else
            assert(isscalar(v), 'Field %s is not scalar', f{j});
            if v == 0, v = 0; end                 % print -0 as 0
            c{j} = sprintf('%.8g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(c, ','));
end
fclose(fid);
end

function writeSettings(fp, C, Gg, VAR, outNames)
% Every setting of the run as name,value text rows (read by alb_analysis.m).
s = {'grid_km', num2str(Gg.agg); 'rows', num2str(Gg.H); 'cols', num2str(Gg.W); ...
     'min_valid_pixels', num2str(Gg.minN); 'years', sprintf('%d-%d', C.years(1), C.years(end)); ...
     'window_days', sprintf('%d-%d', C.window); 'min_frac', num2str(C.minFrac); ...
     'land_frac', num2str(C.landFrac); 'majority_frac', num2str(C.majorityFrac); ...
     'S0_Wm2', num2str(C.S0); 'use_era5', num2str(C.useERA5); 'band_dates', C.bandDates; ...
     'bw_days', num2str(C.opt.bw); 'minwt', num2str(C.opt.minwt); ...
     'pre_window', sprintf('%d-%d', C.opt.pre); 'post_window', sprintf('%d-%d', C.opt.post); ...
     'qhi', num2str(C.opt.qhi); 'qlo', num2str(C.opt.qlo); 'minamp', num2str(C.opt.minamp); ...
     'minobs', num2str(C.opt.minobs); 'layers', strjoin(outNames, ' '); ...
     'albedo_files', fullfile(VAR.alb.dir, VAR.alb.pat); 'ndsi_files', fullfile(VAR.ndsi.dir, VAR.ndsi.pat); ...
     'ssrd_files', fullfile(VAR.ssrd.dir, VAR.ssrd.pat); ...
     'use_alb_qa', num2str(C.useAlbQA); 'min_pct_full', num2str(C.minPctFull); ...
     'albedo_qa_files', fullfile(VAR.albqa.dir, VAR.albqa.pat); ...
     'written', datestr(now, 'yyyy-mm-dd HH:MM'); 'matlab', version};
fid = fopen(fp, 'w');  assert(fid > 0, 'Cannot write %s', fp);
fprintf(fid, 'name,value\n');
for k = 1:size(s, 1), fprintf(fid, '"%s","%s"\n', s{k, 1}, strrep(s{k, 2}, '"', '''')); end
fclose(fid);
end

function checkOutRoot(C, VAR)
% Stop before anything is written if C.outRoot already holds layers prepared
% from other albedo inputs (other files, QA use or C.minPctFull): those results
% are kept rather than silently replaced. Same inputs = a plain re-run.
key = {'albedo_files', fullfile(VAR.alb.dir, VAR.alb.pat); 'use_alb_qa', num2str(C.useAlbQA); ...
       'min_pct_full', num2str(C.minPctFull)};
for a = C.aggs
    fp = fullfile(C.outRoot, sprintf('grid_%dkm', a), 'settings.csv');
    if ~isfile(fp), continue; end
    txt = fileread(fp);
    for k = 1:size(key, 1)
        row = sprintf('"%s","%s"', key{k, 1}, strrep(key{k, 2}, '"', ''''));
        assert(~isempty(strfind(txt, row)), ['%s holds layers prepared from other albedo inputs ' ...
            '(its %s is not %s). Rename that folder (e.g. Prepared_v1) or change C.outRoot.'], ...
            fileparts(fp), key{k, 1}, key{k, 2});                                       %#ok<STREMP>
    end
end
end
