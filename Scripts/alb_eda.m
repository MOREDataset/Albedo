%% ALB_EDA  Exploratory data analysis of the blue-sky Albedo archive of Alaska (1 km, 2003-2025)
%
% Figures (E.figDir, PDF) and tables (E.tabDir, CSV) for the paper's three research questions. Every statistic uses
% the baseline E.base = 2003-2025, the complete years with both Terra and Aqua in MCD43; the other years of E.years
% enter only tab01 and the five-year means (fig07, tab09). Rows of the regional results: the 9 Level-2 ecoregions,
% then All Alaska.
%
% RQ1  How have Alaska's Albedo magnitude, seasonality and spatial distribution changed over 2003-2025?
%   fig01_mean_annual_cycle     daily regional Albedo of every baseline year (7-day smoothed), coloured by its group
%                               of years (E.yearGroups), with each group's mean, the mean of all years and annual
%                               metrics in the legends; fig01_driver_<v>_annual_cycle: the same for NDSI, NDVI, t2m,
%                               SWE and snowfall
%   fig07_five_year_means       five-year Feb-Oct mean Albedo per cell and its difference from the baseline normal
%   fig09_snow_years            statewide daily Albedo and NDSI snow cover, 1 Sep - 30 Aug of the E.snowYears
%   fig16_Albedo_ndsi_annual_cycle   Albedo and NDSI snow cover cycles in the same panels
%   fig17_Albedo_drivers_timeseries  daily Albedo, NDSI, SWE and air temperature, 2003-2025 as one series
%   fig18a-f                    stage 2g: days per year past a threshold and near the summer extreme, with their
%                               Theil-Sen trends: a Albedo below 0.2 (each year's onset, end and duration), b Albedo
%                               at its summer low, c NDSI snow cover below 0.15 (snow-free season), d NDSI at its
%                               summer low, e air temperature above 0 degC (thaw season), f the warmest weeks
%   tab01 dataset by year; tab02 regional monthly climatology; tab09 five-year means by region; tab11 annual Albedo
%   metrics of every region and year, tab11b their trends; tab12 driver annual cycles; tab13 the snow years;
%   tab18 the durations of every region and year, tab18b their trends, tab18c the trends of the onsets and ends
%
% RQ2  How are Albedo anomalies associated with snow, temperature and vegetation across seasons and regions?
%   fig02_driver_anomaly_maps   per cell and season: Spearman rho of the 8-day Albedo and driver anomalies
%   fig03_dependence_rho_xi     Spearman rho and Chatterjee xi of Albedo vs each covariate, by ecoregion
%   tab03 driver anomalies by region, tab03b per cell (GeoTIFF maps above E.maxCsvCells cells), tab03c by
%   permafrost class; tab04 rho and xi statewide, tab04b by region; tab05 Spearman matrices over the pixel-days
%
% RQ3  How sensitive is seasonal Albedo to interannual spring air-temperature variability, and what changes are
%      implied by mid-century CMIP6 warming scenarios? (stage 2h)
%   fig19_temperature_sensitivity  All Alaska: the observed years and, per scenario, the temperature-scaled Albedo
%                               of every year of E.whatIfYears (a box over the models, a line through the yearly means)
%   fig19b_sensitivity_map      b, the Albedo change per degC of April-May air temperature, by ecoregion
%   fig19c_change_maps          the E.whatIfYears Albedo and its change by ecoregion, one map per scenario
%   tab19 fit, checks and estimate of every region and scenario; tab19b every region and year; tab19c the table for
%   the paper; tab19d every region, scenario and year of E.whatIfYears; tab19e every model
%
% METHODS (the comments of the functions named give the details)
%   Cells       k x k pixels of the 1-km grid, k = E.cellKm, with the static layers of Prepared\grid_<k>km; a pixel
%               takes the ecoregion of its cell. Pixel-day thresholds follow from k (cellSettings).
%   Curves      regional daily mean = pooled mean of the region's valid pixels (>= E.minRegDayCov of them valid),
%               centred 7-day mean, 365-day calendar (29 Feb dropped). Annual metrics from each year's curve
%               (cycleMetrics); trends: Theil-Sen slope, Mann-Kendall p, Benjamini-Hochberg q (senKendall, bhFdr).
%   Five-year   per cell and period (E.periods): the mean over Feb-Oct of the period's monthly means.
%   Durations   stage 2g (edaDurations): a season is the run of days around the summer extreme, ended by spells of
%               >= E.durPersist days back across the threshold; the days at peak are counted over the whole year.
%               Trends with the Mann-Kendall test corrected for autocorrelated years (E.durAutocorr).
%   Anomalies   8-day cell means (block k = days 8k-7..8k); anomaly = block mean - the cell's baseline mean of that
%               block. Per cell and season: Spearman rho and Chatterjee xi over the (year, block) pairs
%               (driverMapsFile; >= E.blkMinPairs pairs), mapped as a measure, not tested.
%   Covariates  E.nSample random pixels (seed E.seed): rho and xi of Albedo vs the daily covariates over the
%               pixel-days of a season, and vs the static ones (elevation, slope, northness, latitude, peak NDVI,
%               permafrost class) over the pixels' season mean Albedo; statewide and by ecoregion.
%   RQ3         A = the sunlight-weighted Albedo of the season E.albTempBlocks (18 Feb - 23 Oct; weights: the
%               ERA5-Land SSRD climatology), T = the mean April-May air temperature (E.albTempTBlocks). A = a + b T by
%               least squares over the baseline years, b's E.albTempCI % range from E.albTempNBoot resamples of the
%               years. Year y of E.whatIfYears: A(y) = the mean A + b dT(y), dT(y) = a NEX-GDDP-CMIP6 model's April-May
%               warming of year y (edaAlbTemp). A temperature-scaled estimate, not a projection: it assumes that
%               today's year-to-year relation holds for the forced warming.
%   No Statistics Toolbox: percentiles, ranks, correlations and tests are local functions.
%
% DATA (under E.root; the 1-km EPSG:3338 grid with upper-left corner (-706000, 2375000); the covariate files are one
% row and one column larger and are cropped from the top-left)
%   BlueSky_Daily_TIFFs_v4\BlueSky_Albedo_daily_YYYY.tif        int16 x 1e-4, 1..10000 valid
%   BlueSky_Daily_QA_TIFFs_v4\BlueSky_pctfull_daily_YYYY.tif    % of the value's area from full BRDF inversions
%                                                               (0..100), 255 = nodata
%   SnowCover_Daily_TIFFs\NDSI_SnowCover_MOD10A1_daily_YYYY.tif   0..100 valid
%   NDVI_Daily_TIFFs\MOD09GA_NDVI_daily_YYYY.tif                x 1e-4, -2000..10000 valid
%   ERA5_Daily_TIFFs\ERA5_t2m|swe|snowfall|ssrd_daily_YYYY.tif  ERA5-Land: t2m x 0.01 K (0 = no value), SWE and
%                                                               snowfall x 0.1 mm, SSRD x 0.002 MJ m-2 d-1
%   topography\ALOS_AW3D30_topography.tif   elevation (finite = the AOI), slope and aspect (degrees)
%   F.permafrost   permaice.shp of the Circum-Arctic Map of Permafrost (Brown et al. 2002, NSIDC GGD318)
%   Prepared\grid_<E.cellKm>km\   the static layers of alb_prepare.m
%   E.cmip6Dir     the NEX-GDDP-CMIP6 warming tables of stage 2h (README.md there)
%   Band dates come from the GDAL band descriptions 'dYYYYMMDD'.
%
% RUN
%   Set the USER SETTINGS below and run the script. Stage 1 reads every daily file once and caches per-year sums in
%   E.cacheDir (alb_YYYY.mat, cov_YYYY.mat, sw_YYYY.mat, blk_YYYY_FFF.mat; cached years are skipped unless
%   E.redo = true). Stage 2 analyses the cache and keeps its results in eda_results.mat; stages 2b (eda_extra.mat), 2g
%   and 2h take seconds to minutes. E.figuresOnly = true redraws the figures of E.figures and the tables from the cached
%   results without reading the daily files. Needs the Mapping Toolbox and MATLAB R2021a or later; E.nWorkers > 0 uses
%   the Parallel Computing Toolbox. At 1 km a first full run takes hours (stage 2 alone ~2 h). The history of changes
%   and the test records up to 5 Nov 2026 are in the header of matlab\backup_2026-11-05_before_cleanup\alb_eda.m.
% -------------------------------------------------------------------------

clear; close all; clc;

%% ======================= USER SETTINGS =======================
E.root    = 'C:\Users\aymane.ahajjam\Desktop\Albedo';
E.cellKm  = 1;                % ANALYSIS CELL SIZE (km): k x k pixels of the 1-km grid. It must be a grid
                              % written by alb_prepare.m (C.aggs): the static layers come from
                              % Prepared\grid_<cellKm>km. Every folder, threshold and distance below
                              % follows from it (cellSettings).
E.years   = 2000:2026;        % Albedo years read
E.base    = [2003 2025];      % baseline of every statistic (complete Terra + Aqua years)
E.months  = 2:10;             % Feb..Oct
E.seasons = {'MAM', 60, 151; 'JJA', 152, 243; 'SON', 244, 334};   % name, first and last day of year
E.minCellDaysPx = 4;      % cell mean of a month and year: >= 4 valid days per pixel of the cell
                          % (4 x cellKm^2 pixel-days: 100 at 5 km, 4 at 1 km)
E.minClimYears  = 15;     % cell monthly climatology: >= 15 baseline years
E.minRegShare   = 0.30;   % regional value of a year: >= 30 % of the region's cells observed
E.minRegDayCov  = 0.05;   % regional daily mean: >= 5 % of the region's pixels valid that day
E.minCycleYears = 15;     % mean annual cycle: >= 15 baseline years with a value on that day
E.blkSeasons  = {'All', 1:46; 'MAM', 8:19; 'JJA', 20:30; 'SON', 31:42};   % 8-day blocks of the driver maps
E.blkMinFrac  = 0.1;      % 8-day cell mean: >= 10 % of the cell's pixel-days valid (of 8 x cellKm^2:
                          % 20 of 200 at 5 km, 1 of 8 at 1 km)
E.blkMinYears = 10;       % normal of a cell and block: >= 10 baseline years
E.blkMinPairs = 60;       % driver-anomaly correlation of a cell: >= 60 (year, block) pairs
E.blkChunk    = 5000;     % cells per chunk in the driver maps (memory ~ chunk x 1058 pairs x 8 B x 15, per
                          % worker; the results do not depend on it)
E.blkFileCells = 50000;   % cells per 8-day cache file (stage 2 holds one file's cells for all years:
                          % ~1.3 GB per 50,000 cells)
E.durations   = true;     % stage 2g (fig18a-f, tab18): how many days each year Albedo, NDSI snow cover and air
                          % temperature stay past a threshold and near their summer extreme in every region, and
                          % the trends, from the 7-day smoothed curves of fig01 (false: skip)
E.durAlbedo   = 0.2;      % ... fig18a: the season of Albedo below 0.2, from its spring to its autumn crossing
E.durPeakAlbedo = 0.02;   % ... fig18b, Albedo at its summer low: days below 0.2 and within 0.02 of the year's lowest
E.durNdsi     = 0.15;     % ... fig18c, snow-free season: NDSI snow cover below 0.15, spring to autumn crossing
E.durPeakNdsi = 0.05;     % ... fig18d, NDSI at its summer low: days below 0.15 and within 0.05 of the year's lowest
E.durT2m      = 0;        % ... fig18e, thaw season: 2-m air temperature above 0 degC, spring to autumn crossing
E.durPeakT2m  = 2;        % ... fig18f, the warmest weeks: days within 2 degC of the year's highest air temperature
E.durPersist  = 7;        % ... a season ends at the first spell of >= 7 consecutive days back across its threshold;
                          % shorter spells (a late snowfall, a cold snap, cloud noise) belong to it (1: any day)
E.durMaxGap   = 7;        % ... gaps of <= 7 days are filled linearly; a longer gap leaves the boundary unknown
E.durPivot    = [121 273];   % ... the year's lowest Albedo and NDSI and highest t2m: the extreme of 1 May - 30 Sep
                          % (365-day calendar); a season is the run around its day, the days at peak every day near it
E.durCI       = 95;       % ... intervals of the Theil-Sen slopes (percent)
E.durAutocorr = 'hamedRao';  % ... Mann-Kendall p and slope interval corrected for autocorrelated years (Hamed and Rao
                          % 1998; 'none': the classical test, which assumes independent years)
E.albTemp     = true;     % stage 2h (fig19, fig19b, fig19c, tab19-tab19e; RQ3): the year-to-year sensitivity of the
                          % season's sunlight-weighted Albedo to April-May air temperature, and a temperature-scaled
                          % estimate for every year of E.whatIfYears from the CMIP6 warming of E.cmip6Dir (false: skip)
E.albTempBlocks = 7:37;   % ... the season: 8-day blocks 7-37, 18 Feb - 23 Oct (the longest run of whole blocks with an
                          % Albedo value in every region and year; 92-99 % of the year's sunlight)
E.albTempTBlocks = 12:19; % ... T: mean 2-m air temperature of blocks 12-19, 30 Mar - 1 Jun (April-May, the spring
                          % season of the snow-albedo feedback of Hall & Qu 2006, in whole blocks like the CMIP6 tables)
E.albTempWindow = 2;      % ... check block by block: the Albedo of block k against the temperature of the 2 blocks
                          % (16 days) ending with block k
E.albTempMaxGap = 7;      % ... gaps of <= 7 days in the curves are filled linearly
E.albTempMinYears = 15;   % ... a fit needs >= 15 baseline years with both
E.albTempNBoot = 2000;    % ... resamples of the years (intervals)
E.albTempSeed = 20261027; % ... their seed (rng twister; MATLAB and Octave draw different years from the same seed)
E.albTempCI   = 90;       % ... intervals (percent)
E.cmip6Dir    = fullfile(E.root, 'EDA', 'cmip6');   % ... the CSV tables of the NEX-GDDP-CMIP6 warming
                          % (cmip6_assemble_2031_2050.py; README.md there)
E.whatIfScenarios = {'ssp126', 'ssp245', 'ssp370', 'ssp585'};   % ... scenarios (the four of NEX-GDDP-CMIP6; it has
                          % no SSP1-1.9), two tables each: alb_cmip6_tas_<scenario>_<first>-<last>.csv (block means
                          % of the baseline and of E.whatIfYears) and alb_cmip6_tas_years_<scenario>.csv (every year)
E.whatIfYears = [2031 2050];   % ... the future years [first last]: an estimate for every year and their mean
E.whatIfMinDays = 0.75;   % ... a model is used if every block keeps >= 75 % of its days in both periods (the 360-day
                          % models spread their days over the 365-day year: 7 of 8 days in 3 blocks)
E.periods   = [2000 2004; 2005 2009; 2010 2014; 2015 2019; 2020 2024; 2025 2026];   % fig07 five-year means
E.perMinFrac = 0.5;       % a period's monthly mean needs >= half of its years (3 of 5; 1 of 2 for 2025-26)
E.nSample   = 10000;      % sample pixels of the daily analyses
E.seed      = 20261003;
E.covDays   = [60 334];   % covariates kept for days 60..334 (1 Mar .. 30 Nov)
E.peakMinYears = 5;       % sample-pixel peak NDVI: >= 5 yearly Jul-Aug maxima
E.minPixDays   = 20;      % pixel mean Albedo of a season (static covariates): >= 20 values
E.maxCsvCells  = 100000;  % per-cell table tab03b as CSV up to this many kept cells,
                          % otherwise as GeoTIFF maps on the cell grid (1 km: ~1.5 million cells)
E.outDir   = fullfile(E.root, 'EDA', sprintf('matlab_%dkm', E.cellKm));   % console, cache, figures, tables
E.cacheDir = fullfile(E.outDir, 'cache');
E.figDir   = fullfile(E.outDir, 'figures');
E.tabDir   = fullfile(E.outDir, 'tables');
E.redo        = false;    % true: recompute the cached years
E.figuresOnly = true;    % true: redraw figures and tables from eda_results.mat (and eda_extra.mat when it
                          % matches; no parallel pool)
E.figVisible  = 'on';     % 'off' to draw the figures without showing them
E.nWorkers    = 24;        % parallel workers (Parallel Computing Toolbox) for the per-year passes and the driver
                          % maps; 0 = serial. Without the toolbox the loops run serially. Results do not depend on
                          % it. Each worker needs ~2-3 GB (one 8-day cache file or one year of cells): lower it if
                          % memory runs out
E.metricMinFrac = 0.75;   % annual Albedo metrics (tab11, the fig01 legends): a window of days needs >= 75 %
                          % of its days with a value of the year's 7-day smoothed regional curve
E.regMinN     = 30;       % fig03 by region: rho and xi need >= 30 pixel-days (daily covariates) or pixels (static)
E.snowYears   = [2014 2024];   % optional fig09: statewide daily Albedo and NDSI snow cover from 1 Sep of each of
                          % these years to 30 Aug of the next ([] = no fig09)
E.yearPalette = 'blueRed';  % colours of the fig01 annual cycles: 'blueRed' (green, blue, red) or the light-to-
                          % dark ramps 'tealIndigo', 'viridis', 'batlow' (Crameri), spread evenly over the groups
                          % of E.yearGroups (over the years when it is [])
E.yearGroups  = [2003 2009; 2010 2019; 2020 2025];   % fig01 annual cycles: groups of the baseline years (rows: first
                          % and last year; together every baseline year), one colour each, drawn with their mean
                          % curve and their metrics in the legends; [] = one colour per year and a colour bar
E.groupMinFrac = 0.6;     % a group's mean curve on a day, and its legend metrics, need >= 60 % of its years with a
                          % value (5 of 7, 6 of 10, 4 of 6 years)
E.yearAlpha   = 0.3;      % opacity of the single-year curves when grouped (0-1)
E.bothYears   = false;    % fig16 (Albedo and NDSI in the same panels): true also draws the pale single years
E.bothYLim    = [0 1];    % fig16: limits of the shared axis of Albedo and NDSI snow cover (both 0-1)
E.seriesTLim  = [-40 25]; % fig17 (whole-period series): limits of the air temperature axis (degC, right axis)
E.groupMetrics = {'Albedo', [1 7]; 'ndsi', [1 7]; 'ndvi', 4; 't2m', [1 4]; 'swe', [2 7]; 'snowfall', 1};
                          % annual metrics in the fig01 legends (numbers of cycleMetrics: 1 mean Mar-Oct, 2 spring
                          % peak = maximum Mar-Jun, 4 mean Jun-Aug, 7 half-melt date, 10 snow-free season, ...)
E.cycleLegendFont = 10;   % font size (pt) of the fig01 legends
E.figures     = 19;   % figures drawn: 1 (fig01 and the driver cycles), 2, 3, 7, 9, 16, 17, 18 (fig18a-f) and
                          % 19 (fig19, fig19b, fig19c); 1:19 = all

% input files: folder and file name (%d = year). Paths are joined with fullfile and
% never used as a sprintf format (a backslash there starts an escape sequence).
F.alb      = {fullfile(E.root, 'BlueSky_Daily_TIFFs_v4'),    'BlueSky_Albedo_daily_%d.tif'};
F.qa       = {fullfile(E.root, 'BlueSky_Daily_QA_TIFFs_v4'), 'BlueSky_pctfull_daily_%d.tif'};
F.ndsi     = {fullfile(E.root, 'SnowCover_Daily_TIFFs'),     'NDSI_SnowCover_MOD10A1_daily_%d.tif'};
F.ndvi     = {fullfile(E.root, 'NDVI_Daily_TIFFs'),          'MOD09GA_NDVI_daily_%d.tif'};
F.t2m      = {fullfile(E.root, 'ERA5_Daily_TIFFs'),          'ERA5_t2m_daily_%d.tif'};
F.swe      = {fullfile(E.root, 'ERA5_Daily_TIFFs'),          'ERA5_swe_daily_%d.tif'};
F.snowfall = {fullfile(E.root, 'ERA5_Daily_TIFFs'),          'ERA5_snowfall_daily_%d.tif'};
F.ssrd     = {fullfile(E.root, 'ERA5_Daily_TIFFs'),          'ERA5_ssrd_daily_%d.tif'};
F.topo     = fullfile(E.root, 'topography', 'ALOS_AW3D30_topography.tif');
F.permafrost = 'C:\Users\aymane.ahajjam\Desktop\HighwayStability\data\raw\ggd318_map_circumarctic\permaice.shp';
F.grid     = fullfile(E.root, 'Prepared', sprintf('grid_%dkm', E.cellKm));   % static layers of the cells

%% =============================================================
COV = {'ndsi', 'ndvi', 't2m', 'swe', 'snowfall'};
E = cellSettings(E);                       % thresholds in pixel-days and cells, lags in cells
if E.albTemp                               % stage 2h reads the CMIP6 tables at the end: stop now if one is missing
    assert(numel(E.whatIfYears) == 2, 'E.whatIfYears: [first last] year');   % or lacks a block or a year
    nSc = numel(E.whatIfScenarios);
    for c = [arrayfun(@(s) cmip6Table(E, s, 'period'), 1:nSc, 'UniformOutput', false), ...
            arrayfun(@(s) cmip6Table(E, s, 'years'), 1:nSc, 'UniformOutput', false), {'alb_cmip6_pixels.csv'}]
        assert(isfile(fullfile(E.cmip6Dir, c{1})), ['%s not found: put the CSV tables there ' ...
            '(README.md of E.cmip6Dir, which says how to make them for other E.whatIfYears), or set E.albTemp = ' ...
            'false'], fullfile(E.cmip6Dir, c{1}));
    end
    kNeed = unique([E.albTempTBlocks(:)', E.albTempBlocks(1) - E.albTempWindow + 1 : E.albTempBlocks(end)]);
    yNeed = E.whatIfYears(1):E.whatIfYears(2);
    for s = 1:nSc
        fp = fullfile(E.cmip6Dir, cmip6Table(E, s, 'period'));
        C = readWhatIfCsv(fp, {'model', 'scenario', 'code', 'block', 'tas_base_K', 'tas_future_K', 'n_base_days', ...
            'n_future_days'}, [true true false(1, 6)]);
        assert(all(ismember(kNeed, C.block)), ['stage 2h: %s lacks block(s) %s: put the tables there ' ...
            '(blocks 5-37; README.md of E.cmip6Dir), or set E.albTemp = false'], fp, num2str(setdiff(kNeed, C.block(:)')));
        fp = fullfile(E.cmip6Dir, cmip6Table(E, s, 'years'));
        C = readWhatIfCsv(fp, {'model', 'scenario', 'code', 'block', 'year', 'tas_K', 'n_days'}, [true true false(1, 5)]);
        assert(all(ismember(E.albTempTBlocks, C.block)) && all(ismember(yNeed, C.year)), ['stage 2h: %s lacks ' ...
            'block(s) %s or year(s) %s: put the tables there (README.md of E.cmip6Dir), or set ' ...
            'E.albTemp = false'], fp, num2str(setdiff(E.albTempTBlocks(:)', C.block(:)')), ...
            num2str(setdiff(yNeed, C.year(:)')));
    end
    clear nSc c kNeed yNeed s fp C
end
for d = {E.outDir, E.cacheDir, E.figDir, E.tabDir}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
diary(fullfile(E.outDir, 'alb_eda_console.txt'));  diary on
cleanDiary = onCleanup(@() diary('off'));
fprintf('alb_eda started %s: %d-km cells\n', datestr(now, 'yyyy-mm-dd HH:MM'), E.cellKm);
fprintf('  cell month >= %d valid pixel-days; 8-day cell mean >= %d of %d\n', E.minCellDays, E.blkMinCount, ...
    8 * E.cellKm^2);
resFile = fullfile(E.cacheDir, 'eda_results.mat');

if E.figuresOnly
    nw = 0;                                % figures only: a pool is started below only if a driver cache is missing
else
    nw = poolWorkers(E);                   % parallel workers (0 = serial)
end
if E.figuresOnly
    assert(isfile(resFile), 'E.figuresOnly needs %s from an earlier run', resFile);
    L = load(resFile);  R = L.R;  G = L.G;  clear L
    assert(G.cellKm == E.cellKm, '%s holds %d-km results, E.cellKm is %d', resFile, G.cellKm, E.cellKm);
    smp = samplePixels(G, E);
else
    checkInputs(F, E, COV);
    %% ---------- 0. grids, static layers, sample pixels ----------
    G = staticLayers(F, E);
    [G.pf, G.pfC] = permafrostGrid(F, G, E);   % permafrost extent class of every pixel and cell (cached)
    smp = samplePixels(G, E);
    fprintf('1-km grid %d x %d, AOI %d px; %d-km cells kept %d; sample pixels %d\n', G.H1, G.W1, ...
        nnz(G.aoi), G.cellKm, numel(G.keep), numel(smp.idx));

    %% ---------- 1. per-year pass over the daily files (cached; years in parallel) ----------
    stage1(F, G, smp, E, COV, nw);

    %% ---------- 2. assemble the cache and analyse ----------
    t0 = tic;
    D = assemble(E, G, smp, COV);
    R = edaAnalyze(D, G, E, nw);
    clear D
    save(resFile, 'R', 'G', '-v7');
    fprintf('Analyses done (%.0f s); results kept in %s\n', toc(t0), resFile);
end
printSummary(R, E);

%% ---------- 2b. driver annual cycles, annual Albedo metrics, dependence by region, snow years ----------
% Cheap (minutes) and recomputed on every run, also with E.figuresOnly: the regional daily driver means are
% a cache of their own (drv_YYYY.mat, read once, years in parallel); the rest comes from the cache and R.
% With E.figuresOnly, eda_extra.mat is reused when it was made from the same results and settings.
t0 = tic;
extraFile = fullfile(E.cacheDir, 'eda_extra.mat');  X = [];
if E.figuresOnly && ~E.redo && isfile(extraFile)
    L = load(extraFile);
    if isfield(L.X, 'settings') && isequal(L.X.settings, extraSettings(E, R)), X = L.X; end
    clear L
    if isempty(X), fprintf('  eda_extra.mat was made from other results or settings: recomputing it\n'); end
end
if isempty(X)
    if E.figuresOnly && ~isempty(driverNeed(E)), nw = poolWorkers(E); end
    driverDailyPass(F, G, E, COV, nw);
    X = edaExtra(R, G, E, smp, COV);
    X.settings = extraSettings(E, R);
    save(extraFile, 'X', '-v7');
    fprintf('Driver cycles, annual metrics, dependence by region and snow years done (%.0f s)\n', toc(t0));
else
    fprintf('Driver cycles, annual metrics, dependence by region and snow years: reused %s\n', extraFile);
end
printExtra(X, E);

%% ---------- 2g. how many days each year Albedo, NDSI and air temperature stay past a threshold and near their summer
% extreme, in every region, and the trends ----------
% From the 7-day smoothed regional curves of eda_results.mat (Albedo, R.cycYear) and eda_extra.mat (NDSI, t2m,
% X.drvCycYear); seconds, every run with E.durations (also with E.figuresOnly). DU is kept in
% E.cacheDir\eda_durations.mat for the record
DU = [];
if E.durations
    DU = edaDurations(R, G, E, X);
    save(fullfile(E.cacheDir, 'eda_durations.mat'), 'DU', '-v7');
    printDurations(DU, E);
end

%% ---------- 2h. RQ3: temperature sensitivity of the season's sunlight-weighted Albedo, and temperature-scaled
% estimates for every year of E.whatIfYears (2031-2050) and their mean (the what-if) ----------
% From eda_results.mat (R.cycYear, R.swIn), eda_extra.mat (X.drvCycYear) and the CMIP6 tables in E.cmip6Dir; seconds,
% every run with E.albTemp (also with E.figuresOnly). AT is kept in E.cacheDir\eda_albtemp.mat for the record
AT = [];
if E.albTemp
    AT = edaAlbTemp(R, G, E, X);
    save(fullfile(E.cacheDir, 'eda_albtemp.mat'), 'AT', '-v7');
    printAlbTemp(AT, E);
end

%% ---------- 3. figures and tables ----------
close('all')
drawFigures(R, G, E, yearPalette(E.yearPalette), X);
drawDurationFigures(DU, G, E);
drawAlbTempFigure(AT, G, E);
drawAlbTempMap(AT, G, E);
drawAlbTempChangeMaps(AT, G, E);
writeTables(R, G, E);
writeExtraTables(X, E, G);
writeDurationTables(DU, G, E);
writeAlbTempTables(AT, G, E);
fprintf('alb_eda finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));

%% ======================= LOCAL FUNCTIONS =======================

% ----------------------------- inputs, grids, sample -----------------------------

function fp = fileOf(F, v, yr)
% path of variable v's file of year yr
fp = fullfile(F.(v){1}, sprintf(F.(v){2}, yr));
end

function checkInputs(F, E, COV)
% every daily file must exist before hours of reading
baseYears = E.base(1):E.base(2);
assert(all(ismember(baseYears, E.years)), 'E.base must lie inside E.years');
for v = [{'alb', 'qa', 'ssrd'}, COV]
    yrs = E.years;
    if ismember(v{1}, [COV, {'ssrd'}]), yrs = baseYears; end
    miss = yrs(arrayfun(@(y) ~isfile(fileOf(F, v{1}, y)), yrs));
    assert(isempty(miss), 'Missing %s file(s) for %s in %s', v{1}, mat2str(miss), F.(v{1}){1});
end
assert(isfolder(F.grid), ['Missing %s: run alb_prepare.m with %d in C.aggs, or set E.cellKm to a grid ' ...
    'that Prepared holds'], F.grid, E.cellKm);
for f = {F.topo, F.permafrost, fullfile(F.grid, 'static_ecoregion.tif'), fullfile(F.grid, 'static_land_frac.tif'), ...
         fullfile(F.grid, 'static_lat_deg.tif'), fullfile(F.grid, 'static_elev_mean_m.tif'), ...
         fullfile(F.grid, 'ecoregion_codes.csv'), fullfile(F.grid, 'settings.csv')}
    assert(isfile(f{1}), 'Missing %s', f{1});
end
end

function E = cellSettings(E)
% Settings that depend on the cell size, from the per-pixel and per-area settings: the pixel-day thresholds and the
% classes of the 8-day cache sums and counts.
k = E.cellKm;
assert(isscalar(k) && k >= 1 && k == round(k), 'E.cellKm must be a whole number of km (k x k pixels of the 1-km grid)');
E.minCellDays = ceil(E.minCellDaysPx * k^2 - 1e-9);
E.blkMinCount = ceil(E.blkMinFrac * 8 * k^2 - 1e-9);
assert(E.blkMinCount >= 1, 'E.blkMinFrac gives an 8-day threshold below one pixel-day');
% 8-day cache: sums of whole raw values (|raw| <= 32767) over <= 8 k^2 pixel-days; single is exact
% below 2^24, uint8 holds counts up to 255
if 8 * k^2 * 32767 < 2^24, E.blkSumClass = 'single'; else, E.blkSumClass = 'double'; end
if 8 * k^2 <= 255, E.blkCntClass = 'uint8'; else, E.blkCntClass = 'uint16'; end
E.blkName = 'blk';        % 8-day cache files blk_YYYY_FFF.mat
end

function st = readGridSettings(fp)
% settings.csv of alb_prepare.m ("name","value" rows) -> struct of strings
txt = strtrim(fileread(fp));
lines = regexp(txt, '\r?\n', 'split');
st = struct();
for k = 2:numel(lines)
    f = regexp(strtrim(lines{k}), '^"?([^",]*)"?,"?(.*?)"?$', 'tokens', 'once');
    if ~isempty(f), st.(f{1}) = f{2}; end
end
assert(isfield(st, 'grid_km') && isfield(st, 'land_frac'), '%s: expected the rows grid_km and land_frac', fp);
end

function G = staticLayers(F, E)
% 1-km grid, AOI, elevation; cells of E.cellKm (G.agg x G.agg pixels from the top-left);
% region code of every 1-km pixel (1..9 = ecoregion of its cell, 10 = AOI in a cell
% without one, 0 = outside the AOI); kept cells (an ecoregion or >= land_frac of
% settings.csv land; column-major order) and their region.
G.origin = [-706000 2375000];  G.cell = 1000;  G.cellKm = E.cellKm;  G.agg = E.cellKm * 1000 / G.cell;
G.gridDir = F.grid;
ref = fileOf(F, 'alb', E.years(1));
gi = georasterinfo(ref);
G.H1 = gi.RasterSize(1);  G.W1 = gi.RasterSize(2);
checkOrigin(ref, G.origin, G.cell);
G.HC = floor(G.H1 / G.agg);  G.WC = floor(G.W1 / G.agg);
st = readGridSettings(fullfile(F.grid, 'settings.csv'));
assert(str2double(st.grid_km) == E.cellKm, '%s is a %s-km grid, E.cellKm is %d', F.grid, st.grid_km, E.cellKm);
G.landFrac = str2double(st.land_frac);
[ecoC, RC] = readgeoraster(fullfile(F.grid, 'static_ecoregion.tif'), 'OutputType', 'double');
lfC  = readgeoraster(fullfile(F.grid, 'static_land_frac.tif'), 'OutputType', 'double');
latC = readgeoraster(fullfile(F.grid, 'static_lat_deg.tif'), 'OutputType', 'double');
elC  = readgeoraster(fullfile(F.grid, 'static_elev_mean_m.tif'), 'OutputType', 'double');
assert(isequal(size(ecoC), [G.HC G.WC]) && isequal(size(lfC), [G.HC G.WC]) && isequal(size(latC), [G.HC G.WC]) && ...
    isequal(size(elC), [G.HC G.WC]), ...
    'The layers in %s are not %d x %d cells', F.grid, G.HC, G.WC);
assert(RC.XWorldLimits(1) == G.origin(1) && RC.YWorldLimits(2) == G.origin(2) && ...
    RC.CellExtentInWorldX == G.agg * G.cell, 'The layers in %s are not on the %d-km grid of the Albedo files', ...
    F.grid, E.cellKm);
ecoC(~isfinite(ecoC)) = 0;
codes = readCodes(fullfile(F.grid, 'ecoregion_codes.csv'));
assert(isequal(codes.code(:)', 1:9), 'Expected ecoregion codes 1..9 in ecoregion_codes.csv');
G.regNames = codes.level2;  G.regL1 = codes.level1;
G.order = [3 5 4 8 1 9 6 7 2];                         % plots: Polar, Boreal, Maritime
G.short = {'Alaska Range Trans.', 'Aleutian Meadows', 'Arctic Tundra', 'Bering Taiga', 'Bering Tundra', ...
           'Coast Mtns Trans.', 'Coastal Rainforests', 'Intermontane Boreal', 'Pacific Mtns Trans.'};
[topo, Rt] = readgeoraster(F.topo, 'OutputType', 'single');
assert(Rt.XWorldLimits(1) == G.origin(1) && Rt.YWorldLimits(2) == G.origin(2) && Rt.CellExtentInWorldX == G.cell, ...
    'The topography file is not on the 1-km grid of the Albedo files');
assert(size(topo, 1) >= G.H1 && size(topo, 2) >= G.W1, 'Unexpected size of the topography file');
assert(size(topo, 3) >= 3, ['The topography file must hold elevation, slope and aspect (3 bands); ' ...
    're-export it as described in DATA']);
G.elev   = double(topo(1:G.H1, 1:G.W1, 1));               % m
G.slope  = double(topo(1:G.H1, 1:G.W1, 2));               % degrees
G.aspect = double(topo(1:G.H1, 1:G.W1, 3));               % degrees clockwise from north
clear topo
G.aoi = isfinite(G.elev);
code = zeros(G.H1, G.W1);
code(1:G.HC*G.agg, 1:G.WC*G.agg) = kron(ecoC, ones(G.agg));
code(G.aoi & code == 0) = 10;
code(~G.aoi) = 0;
G.code = int8(code);
G.ecoC = ecoC;  G.lfC = lfC;  G.latC = latC;  G.elevC = elC;
G.keep = find(ecoC(:) > 0 | lfC(:) >= G.landFrac - 1e-9);
G.regC = ecoC(G.keep);
G.npix = arrayfun(@(k) nnz(code == k), 0:10);          % pixels per code 0..10
grpNames = {'Polar', 'Boreal', 'Maritime'};
G.grpNames = grpNames;
G.group = cellfun(@(s) find(strcmp(grpNames, s)), G.regL1);   % Level-1 group 1..3 of region code 1..9
G.x1 = (G.origin(1) + ((1:G.W1) - 0.5) * G.cell) / 1e3;          % km, pixel centres
G.y1 = (G.origin(2) - ((1:G.H1) - 0.5) * G.cell) / 1e3;
G.xc = (G.origin(1) + ((1:G.WC) - 0.5) * G.agg * G.cell) / 1e3;  % km, cell centres
G.yc = (G.origin(2) - ((1:G.HC) - 0.5) * G.agg * G.cell) / 1e3;
G.cities = {'Utqiagvik', -102402, 2367965; 'Nome', -545163, 1662234; 'Fairbanks', 297506, 1667290; ...
            'Bethel', -420091, 1225219; 'Anchorage', 219352, 1255515; 'Juneau', 1127961, 1089735};
end

function c = readCodes(fp)
% ecoregion_codes.csv of alb_prepare.m: code, level2, level1, ... (text fields quoted)
txt = strtrim(fileread(fp));
lines = regexp(txt, '\r?\n', 'split');
hdr = strrep(strsplit(strtrim(lines{1}), ','), '"', '');
iC = find(strcmp(hdr, 'code'));  i2 = find(strcmp(hdr, 'level2'));  i1 = find(strcmp(hdr, 'level1'));
assert(~isempty(iC) && ~isempty(i2) && ~isempty(i1), '%s: expected columns code, level2, level1', fp);
n = numel(lines) - 1;
c.code = zeros(n, 1);  c.level2 = cell(1, n);  c.level1 = cell(1, n);
for k = 1:n
    f = regexp(strtrim(lines{k + 1}), '("[^"]*"|[^,]*)', 'match');
    f = strrep(f, '"', '');
    c.code(k) = str2double(f{iC});  c.level2{k} = f{i2};  c.level1{k} = f{i1};
end
[c.code, o] = sort(c.code);  c.level2 = c.level2(o);  c.level1 = c.level1(o);
end

function smp = samplePixels(G, E)
% E.nSample random AOI pixels (column-major linear indices into the 1-km grid), drawn
% once with seed E.seed and kept in the cache, so every year uses the same pixels.
fs = fullfile(E.cacheDir, 'samples.mat');
if isfile(fs)
    smp = load(fs);
    ok = smp.n == E.nSample && smp.seed == E.seed && isequal(smp.gridSize, [G.H1 G.W1]) && all(G.aoi(smp.idx));
    assert(ok, ['%s was drawn with other settings or for another grid: delete it and the ' ...
        'alb_*.mat and cov_*.mat files of the cache, then run again'], fs);
    return
end
pool = find(G.aoi);
assert(numel(pool) >= E.nSample, 'The AOI holds fewer than E.nSample pixels');
rng(E.seed, 'twister');
idx = sort(pool(randperm(numel(pool), E.nSample)));
smp = struct('idx', idx(:), 'n', E.nSample, 'seed', E.seed, 'gridSize', [G.H1 G.W1]);
save(fs, '-struct', 'smp');
end

function ok = cacheOk(fp, need, cellKm)
% true if the cache file exists, holds the variables this version writes and was made for
% cells of cellKm
ok = false;
if ~isfile(fp), return; end
w = whos('-file', fp);
ok = all(ismember(need, {w.name}));
if ok, ok = getfield(load(fp, 'cellKm'), 'cellKm') == cellKm; end                  %#ok<GFLD>
end

function [s, loc] = tileSamples(smp, G, rr, cc)
% sample numbers inside a tile and their pixel index within the tile (column-major)
[r, c] = ind2sub([G.H1 G.W1], smp.idx);
s = find(r >= rr(1) & r <= rr(end) & c >= cc(1) & c <= cc(end));
loc = (c(s) - cc(1)) * numel(rr) + (r(s) - rr(1)) + 1;
end

% ----------------------------- stage 1: per-year passes -----------------------------

function S = passAlbedo(yr, F, G, smp, E)
% One year of Albedo + QA -> (months E.months; region codes 0..10 = rows 1..11;
% days = actual day of year 1..366):
%   msum / mcnt [H1 x W1 x nMonths]  per 1-km pixel and month: sum of the valid
%                                    Albedo values, number of valid days
%   rsum / rcnt / rmag [11 x 366]    per region code and day: sum of Albedo, valid
%                                    values, sum of the magnitude-inversion share
%   samp [nSample x 366]             raw daily Albedo at the sample pixels (-32768 = none)
fa = fileOf(F, 'alb', yr);  fq = fileOf(F, 'qa', yr);
checkOrigin(fa, G.origin, G.cell);  checkOrigin(fq, G.origin, G.cell);
dA = tagDates(fa);  dQ = tagDates(fq);
[dn, ia, iq] = intersect(dA, dQ);                      % band dates both files hold
assert(~isempty(dn), '%s and %s share no band date', fa, fq);
if numel(dn) < max(numel(dA), numel(dQ))               % e.g. a partial year exported on different days
    fprintf('  note: %s (%d bands) and %s (%d bands) share %d band dates; only those are used\n', ...
        fa, numel(dA), fq, numel(dQ), numel(dn));
end
dv = datevec(dn);
assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fa, yr);
mon = dv(:, 2);  doy = round(dn - datenum(yr, 1, 0));
nB = numel(dn);  nM = numel(E.months);
Mon = double(mon == E.months(:)');                     % nB x nM, band -> month
S = struct('msum', zeros(G.H1, G.W1, nM, 'single'), 'mcnt', zeros(G.H1, G.W1, nM, 'uint8'), ...
    'rsum', zeros(11, 366), 'rcnt', zeros(11, 366), 'rmag', zeros(11, 366), ...
    'samp', repmat(int16(-32768), numel(smp.idx), 366), 'months', E.months, ...
    'nB', nB, 'nA', numel(dA), 'nQ', numel(dQ), 'valid', 0, 'mismatch', 0, 'first', dn(1), 'last', dn(end), ...
    'cellKm', G.cellKm, 'cacheVersion', 3);
Fa = openTiff(fa);  closeA = onCleanup(@() closeTiff(Fa));
Fq = openTiff(fq);  closeQ = onCleanup(@() closeTiff(Fq));
assert(Fq.tw == Fa.tw && Fq.th == Fa.th, '%s and %s are tiled differently', fa, fq);
assert(Fa.H == G.H1 && Fa.W == G.W1 && Fq.H == G.H1 && Fq.W == G.W1, '%s or %s is not %d x %d', fa, fq, G.H1, G.W1);
for ty = 0:ceil(G.H1 / Fa.th) - 1
    rr = ty*Fa.th + 1 : min((ty+1)*Fa.th, G.H1);
    for tx = 0:ceil(G.W1 / Fa.tw) - 1
        cc = tx*Fa.tw + 1 : min((tx+1)*Fa.tw, G.W1);
        A = readTile(Fa, rr, cc);  A = A(:, ia);
        Q = readTile(Fq, rr, cc);  Q = Q(:, iq);
        T = tileAlbedo(A, Q, double(reshape(G.code(rr, cc), [], 1)), Mon);
        S.rsum(:, doy) = S.rsum(:, doy) + T.rsum;
        S.rcnt(:, doy) = S.rcnt(:, doy) + T.rcnt;
        S.rmag(:, doy) = S.rmag(:, doy) + T.rmag;
        S.msum(rr, cc, :) = reshape(T.msum, numel(rr), numel(cc), nM);
        S.mcnt(rr, cc, :) = uint8(reshape(T.mcnt, numel(rr), numel(cc), nM));
        S.valid = S.valid + T.valid;  S.mismatch = S.mismatch + T.mismatch;
        [s, loc] = tileSamples(smp, G, rr, cc);
        if ~isempty(s)
            v = A(loc, :);  v(~(v >= 1 & v <= 10000)) = -32768;
            S.samp(s, doy) = int16(v);
        end
    end
end
assert(S.mismatch <= 1e-3 * max(S.valid, 1), ['%s: %d valid Albedo values have no QA value - ' ...
    'is %s the companion of this file?'], fa, S.mismatch, fq);
end

function T = tileAlbedo(A, Q, cd, Mon)
% Sums of one tile: A Albedo and Q QA (pixels x bands), cd region code 0..10 of
% every pixel, Mon band -> month indicator (bands x nM).
m = A >= 1 & A <= 10000;
qok = Q <= 100;
T.valid = nnz(m);  T.mismatch = nnz(m & ~qok);
X = single(A) * single(1e-4);  X(~m) = 0;
Mf = single(m);
Gm = 1 - single(Q) / 100;  Gm(~(m & qok)) = 0;
P = sparse(cd + 1, 1:numel(cd), 1, 11, numel(cd));    % region code -> row
T.rsum = full(P * double(X));  T.rcnt = full(P * double(Mf));  T.rmag = full(P * double(Gm));
T.msum = X * single(Mon);                              % pixels x nM
T.mcnt = Mf * single(Mon);
end

function S = passCovariates(yr, F, G, smp, COV)
% raw daily values of the covariates at the sample pixels [nSample x 366], -32768 = none
for v = COV
    fp = fileOf(F, v{1}, yr);
    checkOrigin(fp, G.origin, G.cell);
    dn = tagDates(fp);
    dv = datevec(dn);
    assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fp, yr);
    doy = round(dn - datenum(yr, 1, 0));
    out = repmat(int16(-32768), numel(smp.idx), 366);
    Fk = openTiff(fp);  ck = onCleanup(@() closeTiff(Fk));
    assert(Fk.H >= G.H1 && Fk.W >= G.W1, '%s is smaller than the Albedo grid', fp);
    for ty = 0:ceil(G.H1 / Fk.th) - 1
        rr = ty*Fk.th + 1 : min((ty+1)*Fk.th, G.H1);
        for tx = 0:ceil(G.W1 / Fk.tw) - 1
            cc = tx*Fk.tw + 1 : min((tx+1)*Fk.tw, G.W1);
            [s, loc] = tileSamples(smp, G, rr, cc);
            if isempty(s), continue; end
            a = readTile(Fk, rr, cc);
            out(s, doy) = int16(a(loc, :));
        end
    end
    S.(v{1}) = out;
    clear ck                                           % closes the file
end
S.cellKm = G.cellKm;  S.cacheVersion = 3;
end

function S = passShortwave(yr, F, G)
% ERA5-Land SSRD of one year per region code (rows 1..11 = codes 0..10) and actual day of year (1..366), over the
% 1-km pixels: inAll = the sum of SSRD (raw) over the ERA5-Land pixels with a valid value, nAll = their number, on
% the days that the year's Albedo file also holds. SSRD raw x 0.002 = MJ m-2 d-1, valid 0..25000. ERA5-Land has no
% value over sea: such pixels hold -32768 on most days but 0 on a few (an export artefact), so a pixel counts as
% ERA5-Land only if it has a value on every day of the year (zeroFilled counts the AOI pixel-days set to 0 off land).
fa = fileOf(F, 'alb', yr);  fs = fileOf(F, 'ssrd', yr);
checkOrigin(fs, G.origin, G.cell);
dA = tagDates(fa);  dS = tagDates(fs);
[dn, ~, is] = intersect(dA, dS);
assert(~isempty(dn), '%s and %s share no band date', fa, fs);
dv = datevec(dn);
assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fs, yr);
doy = round(dn - datenum(yr, 1, 0));
z = zeros(11, 366);
S = struct('inAll', z, 'nAll', z, 'zeroFilled', 0, 'nLand', 0, 'nDays', numel(dn), 'nA', numel(dA), 'nS', numel(dS), ...
    'cellKm', G.cellKm, 'cacheVersion', 1);
Fs = openTiff(fs);  closeS = onCleanup(@() closeTiff(Fs));
assert(Fs.H >= G.H1 && Fs.W >= G.W1, '%s is smaller than the Albedo grid', fs);
for ty = 0:ceil(G.H1 / Fs.th) - 1
    rr = ty*Fs.th + 1 : min((ty+1)*Fs.th, G.H1);
    for tx = 0:ceil(G.W1 / Fs.tw) - 1
        cc = tx*Fs.tw + 1 : min((tx+1)*Fs.tw, G.W1);
        cd = double(reshape(G.code(rr, cc), [], 1));
        if ~any(cd), continue; end
        Sall = readTile(Fs, rr, cc);
        land = all(Sall ~= -32768, 2) & cd > 0;
        Sv = double(Sall(:, is));  clear Sall
        S.zeroFilled = S.zeroFilled + nnz(Sv == 0 & ~land & cd > 0);
        S.nLand = S.nLand + nnz(land);
        ms = land & Sv >= 0 & Sv <= 25000;
        P = sparse(cd + 1, 1:numel(cd), 1, 11, numel(cd));          % region code -> row
        x = Sv;  x(~ms) = 0;  S.inAll(:, doy) = S.inAll(:, doy) + full(P * x);
        S.nAll(:, doy) = S.nAll(:, doy) + full(P * double(ms));
    end
end
end

function S = passBlocks(yr, F, G, E)
% Kept cells x the 46 8-day blocks of the year [nK x 46 x 6]: sum of the valid raw values of
% Albedo, NDSI, NDVI, t2m, SWE and snowfall over the cell's AOI pixel-days (s, E.blkSumClass)
% and their number (n, E.blkCntClass). Validity as in the analyses: Albedo 1..10000, NDSI
% 0..100, NDVI -2000..10000, t2m not -32768 or 0 (no ERA5-Land value), SWE and snowfall only
% where t2m is valid, SWE 32767 (saturated) excluded. The sums are of whole raw values and
% exact in the class chosen by cellSettings. Each tile adds to the cells it covers only.
blocks = 1:46;  nb = numel(blocks);  nK = numel(G.keep);
vars = {'alb', 'ndsi', 'ndvi', 't2m', 'swe', 'snowfall'};
S = struct('s', zeros(nK, nb, 6, E.blkSumClass), 'n', zeros(nK, nb, 6, E.blkCntClass), 'blocks', blocks, ...
    'vars', {vars}, 'cellKm', G.cellKm, 'cacheVersion', 4);
kix = keptIndex(G);                                                % pixel -> kept cell (0 = none)
for v = 1:3                                                        % Albedo, NDSI, NDVI
    fp = fileOf(F, vars{v}, yr);
    checkOrigin(fp, G.origin, G.cell);
    [Ind, use] = blockIndicator(tagDates(fp), blocks);
    Fk = openTiff(fp);  ck = onCleanup(@() closeTiff(Fk));
    assert(Fk.H >= G.H1 && Fk.W >= G.W1, '%s is smaller than the Albedo grid', fp);
    for ty = 0:ceil(G.H1 / Fk.th) - 1
        rr = ty*Fk.th + 1 : min((ty+1)*Fk.th, G.H1);
        for tx = 0:ceil(G.W1 / Fk.tw) - 1
            cc = tx*Fk.tw + 1 : min((tx+1)*Fk.tw, G.W1);
            k = kix(rr, cc);  k = k(:);
            if ~any(k), continue; end
            a = double(readTile(Fk, rr, cc));  a = a(:, use);
            switch v
                case 1, m = a >= 1 & a <= 10000;
                case 2, m = a >= 0 & a <= 100;
                case 3, m = a >= -2000 & a <= 10000;
            end
            a(~m) = 0;
            [u, ds, dc] = tileBlocks(k, a, m, Ind);
            S.s(u, :, v) = S.s(u, :, v) + cast(ds, E.blkSumClass);
            S.n(u, :, v) = S.n(u, :, v) + cast(dc, E.blkCntClass);
        end
    end
    clear ck                                                       % closes the file
end
fe = cellfun(@(x) fileOf(F, x, yr), vars(4:6), 'UniformOutput', false);    % ERA5: t2m, SWE, snowfall
dn = tagDates(fe{1});
for q = 1:3
    checkOrigin(fe{q}, G.origin, G.cell);
    if q > 1, assert(isequal(tagDates(fe{q}), dn), '%s and %s hold different dates', fe{1}, fe{q}); end
end
[Ind, use] = blockIndicator(dn, blocks);
Ft = openTiff(fe{1});  c1 = onCleanup(@() closeTiff(Ft));
Fw = openTiff(fe{2});  c2 = onCleanup(@() closeTiff(Fw));
Ff = openTiff(fe{3});  c3 = onCleanup(@() closeTiff(Ff));
assert(Fw.th == Ft.th && Ff.th == Ft.th && Fw.tw == Ft.tw && Ff.tw == Ft.tw, 'The ERA5 files are tiled differently');
for ty = 0:ceil(G.H1 / Ft.th) - 1
    rr = ty*Ft.th + 1 : min((ty+1)*Ft.th, G.H1);
    for tx = 0:ceil(G.W1 / Ft.tw) - 1
        cc = tx*Ft.tw + 1 : min((tx+1)*Ft.tw, G.W1);
        k = kix(rr, cc);  k = k(:);
        if ~any(k), continue; end
        t = double(readTile(Ft, rr, cc));  t = t(:, use);
        w = double(readTile(Fw, rr, cc));  w = w(:, use);
        f = double(readTile(Ff, rr, cc));  f = f(:, use);
        tok = t ~= -32768 & t ~= 0;
        mw = tok & w ~= -32768 & w ~= 32767;  mf = tok & f ~= -32768;
        t(~tok) = 0;  w(~mw) = 0;  f(~mf) = 0;
        X = {t, w, f};  M = {tok, mw, mf};
        for q = 1:3
            [u, ds, dc] = tileBlocks(k, X{q}, M{q}, Ind);
            S.s(u, :, 3 + q) = S.s(u, :, 3 + q) + cast(ds, E.blkSumClass);
            S.n(u, :, 3 + q) = S.n(u, :, 3 + q) + cast(dc, E.blkCntClass);
        end
    end
end
end

function kix = keptIndex(G)
% kept-cell number (1..numel(G.keep)) of every 1-km AOI pixel inside a kept cell, 0 elsewhere [H1 x W1]
nK = numel(G.keep);
pos = zeros(G.HC * G.WC, 1);  pos(G.keep) = 1:nK;                  % cell -> kept number
[r1, c1] = ndgrid(1:G.H1, 1:G.W1);
inC = r1 <= G.HC * G.agg & c1 <= G.WC * G.agg & G.aoi;
kix = zeros(G.H1, G.W1);
kix(inC) = pos(sub2ind([G.HC G.WC], ceil(r1(inC) / G.agg), ceil(c1(inC) / G.agg)));
end

function [Ind, use] = blockIndicator(dn, blocks)
% band -> 8-day block indicator (bands used x numel(blocks)); block k of the 365-day
% calendar = days 8k-7..8k, 29 Feb in no block
dv = datevec(dn);  yr = dv(:, 1);
doy = round(dn(:) - datenum(yr, 1, 0));
leap = (mod(yr, 4) == 0 & mod(yr, 100) ~= 0) | mod(yr, 400) == 0;
b = floor((doy - (leap & doy > 59) - 1) / 8) + 1;
b(leap & doy == 60) = 0;
Ind = double(b == blocks(:)');
use = find(any(Ind, 2));
Ind = Ind(use, :);
end

function [u, ds, dc] = tileBlocks(k, a, m, Ind)
% one tile's block sums per kept cell: a (pixels x bands, invalid set to 0), m (valid), k the
% kept cell of every pixel (0 = none) -> the cells u it covers, their sums ds and counts dc
% (whole numbers in double, exact)
in = k > 0;
[u, ~, j] = unique(k(in));
P = sparse(j, 1:numel(j), 1, numel(u), numel(j));
ds = full(P * (a(in, :) * Ind));
dc = full(P * (double(m(in, :)) * Ind));
end

function fp = blockFile(E, yr, f)
% 8-day cache file f of year yr (cells (f-1) * E.blkFileCells + 1 ...): E.blkName 'blk' (every valid Albedo value
% and the drivers, passBlocks)
fp = fullfile(E.cacheDir, sprintf('%s_%d_%03d.mat', E.blkName, yr, f));
end

function saveBlocks(E, G, yr, S)
% the year's 8-day sums and counts in files of E.blkFileCells kept cells
nK = numel(G.keep);
for f = 1:ceil(nK / E.blkFileCells)
    rows = (f - 1) * E.blkFileCells + 1 : min(f * E.blkFileCells, nK);
    C = struct('s', S.s(rows, :, :), 'n', S.n(rows, :, :), 'rows', [rows(1) rows(end)], 'nK', nK, ...
        'blocks', S.blocks, 'vars', {S.vars}, 'cellKm', S.cellKm, 'cacheVersion', S.cacheVersion);
    save(blockFile(E, yr, f), '-struct', 'C', '-v7');
end
end

function ok = blockCacheOk(E, G, yr)
% true if every 8-day cache file of the year exists and matches the cells and their split
nK = numel(G.keep);  ok = true;
for f = 1:ceil(nK / E.blkFileCells)
    fp = blockFile(E, yr, f);
    if ~cacheOk(fp, {'s', 'n', 'rows', 'nK', 'blocks', 'vars', 'cellKm', 'cacheVersion'}, E.cellKm)
        ok = false;  return
    end
    L = load(fp, 'rows', 'nK', 'blocks', 'cacheVersion');
    if ~(isequal(L.rows, [(f - 1) * E.blkFileCells + 1, min(f * E.blkFileCells, nK)]) && L.nK == nK && ...
            isequal(L.blocks, 1:46) && L.cacheVersion == 4)
        ok = false;  return
    end
end
end

% ----------------------------- parallel stage 1 and the regional daily driver means -----------------------------

function nw = poolWorkers(E)
% number of parallel workers for the parfor loops: E.nWorkers if the Parallel Computing Toolbox is there
% (a local process pool of that size is started or reused), otherwise 0 (parfor runs as a for loop)
nw = 0;
if E.nWorkers <= 0, return; end
if ~(license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')))
    fprintf('  no Parallel Computing Toolbox: the loops run serially\n');
    return
end
p = gcp('nocreate');
if ~isempty(p) && p.NumWorkers ~= E.nWorkers, delete(p);  p = []; end
if isempty(p)
    try
        p = parpool('Processes', E.nWorkers);               % R2022b+
    catch
        p = parpool('local', E.nWorkers);
    end
end
nw = p.NumWorkers;
fprintf('  parallel pool: %d workers\n', nw);
end

function stage1(F, G, smp, E, COV, nw)
% the per-year passes over the daily files, one year per worker; each year writes its own cache files.
% If the parallel loop fails (e.g. a worker cannot reach a function), it is run again serially: the years
% already cached are skipped, and a real error repeats in the serial run.
yrs = E.years;
try
    parfor (iy = 1:numel(yrs), nw)
        stage1Year(yrs(iy), F, G, smp, E, COV);
    end
catch err
    if nw == 0, rethrow(err); end
    fprintf('  the parallel per-year pass failed (%s): running it serially\n', err.message);
    for iy = 1:numel(yrs), stage1Year(yrs(iy), F, G, smp, E, COV); end
end
end

function stage1Year(yr, F, G, smp, E, COV)
% the cache files of one year
baseYears = E.base(1):E.base(2);
fa = fullfile(E.cacheDir, sprintf('alb_%d.mat', yr));
if E.redo || ~cacheOk(fa, {'msum', 'mcnt', 'rsum', 'samp', 'months', 'cellKm', 'cacheVersion'}, E.cellKm)
    if isfile(fa), fprintf('%d: replacing the cache file of an earlier alb_eda version\n', yr); end
    t0 = tic;
    S = passAlbedo(yr, F, G, smp, E);
    saveStruct(fa, S);
    fprintf('%d Albedo + QA: %d bands (%s to %s), %d valid values, %d without a QA value (%.0f s)\n', ...
        yr, S.nB, datestr(S.first, 'yyyy-mm-dd'), datestr(S.last, 'yyyy-mm-dd'), S.valid, S.mismatch, toc(t0));
end
if ~ismember(yr, baseYears), return; end
fc = fullfile(E.cacheDir, sprintf('cov_%d.mat', yr));
if E.redo || ~cacheOk(fc, [COV, {'cellKm', 'cacheVersion'}], E.cellKm)
    t0 = tic;
    S = passCovariates(yr, F, G, smp, COV);
    saveStruct(fc, S);
    fprintf('%d covariates at the sample pixels (%.0f s)\n', yr, toc(t0));
end
fw = fullfile(E.cacheDir, sprintf('sw_%d.mat', yr));
if E.redo || ~cacheOk(fw, {'inAll', 'nAll', 'cellKm', 'cacheVersion'}, E.cellKm)
    t0 = tic;
    S = passShortwave(yr, F, G);
    saveStruct(fw, S);
    fprintf('%d shortwave: %d days, %d ERA5-Land pixels, %d off-land pixel-days set to 0 ignored (%.0f s)\n', ...
        yr, S.nDays, S.nLand, S.zeroFilled, toc(t0));
end
if E.redo || ~blockCacheOk(E, G, yr)
    t0 = tic;
    S = passBlocks(yr, F, G, E);
    saveBlocks(E, G, yr, S);
    fprintf('%d 8-day cell means of Albedo and drivers, %d cache files (%.0f s)\n', yr, ...
        ceil(numel(G.keep) / E.blkFileCells), toc(t0));
end
end

function saveStruct(fp, S)
% save the fields of S as variables (callable inside parfor, where save cannot be used directly)
save(fp, '-struct', 'S', '-v7');
end

function yrs = driverYears(E)
% years of the regional daily driver means: the baseline (driver cycles) and the years of the snow-year
% windows (E.snowYears: 1 Sep of year y to 30 Aug of y + 1)
yrs = E.base(1):E.base(2);
if ~isempty(E.snowYears)
    sy = E.snowYears(:)';
    yrs = unique([yrs, sy, sy + 1]);
    assert(all(ismember([sy, sy + 1], E.years)), 'E.snowYears needs the Albedo of years %s, outside E.years', ...
        mat2str(setdiff([sy, sy + 1], E.years)));
end
end

function need = driverNeed(E)
% years of driverYears whose drv_YYYY.mat is missing or out of date (all of them with E.redo)
yrs = driverYears(E);
need = yrs(arrayfun(@(y) E.redo || ~cacheOk(drvFile(E, y), {'sum', 'cnt', 'vars', 'cellKm', 'cacheVersion'}, E.cellKm), yrs));
end

function driverDailyPass(F, G, E, COV, nw)
% drv_YYYY.mat for every year of driverYears that is not cached yet (years in parallel)
need = driverNeed(E);
if isempty(need), return; end
for v = COV                                                         % every file before a long read
    miss = need(arrayfun(@(y) ~isfile(fileOf(F, v{1}, y)), need));
    assert(isempty(miss), 'Missing %s file(s) for %s in %s', v{1}, mat2str(miss), F.(v{1}){1});
end
try
    parfor (iy = 1:numel(need), nw)
        driverDailyYear(need(iy), F, G, E, COV);
    end
catch err
    if nw == 0, rethrow(err); end
    fprintf('  the parallel driver pass failed (%s): running it serially\n', err.message);
    for iy = 1:numel(need), driverDailyYear(need(iy), F, G, E, COV); end
end
end

function fp = drvFile(E, yr)
fp = fullfile(E.cacheDir, sprintf('drv_%d.mat', yr));
end

function driverDailyYear(yr, F, G, E, COV)
fp = drvFile(E, yr);
if ~E.redo && cacheOk(fp, {'sum', 'cnt', 'vars', 'cellKm', 'cacheVersion'}, E.cellKm), return; end
t0 = tic;
S = passDriverDaily(yr, F, G, COV);
saveStruct(fp, S);
fprintf('%d regional daily driver means (%.0f s)\n', yr, toc(t0));
end

function S = passDriverDaily(yr, F, G, COV)
% Regional daily sums of the drivers COV (NDSI, NDVI, t2m, SWE, snowfall) over the 1-km pixels, as rsum /
% rcnt of the Albedo: sum [11 x 366 x 5] of the valid raw values per region code 0..10 (rows 1..11) and actual
% day of year, cnt the number of valid pixel-days. Validity as in passBlocks: NDSI 0..100, NDVI -2000..10000,
% t2m not -32768 or 0 (no ERA5-Land value), SWE and snowfall only where t2m is valid, SWE 32767 (saturated)
% excluded.
assert(isequal(COV, {'ndsi', 'ndvi', 't2m', 'swe', 'snowfall'}), 'passDriverDaily expects COV = ndsi, ndvi, t2m, swe, snowfall');
S = struct('sum', zeros(11, 366, 5), 'cnt', zeros(11, 366, 5), 'vars', {COV}, 'nDays', zeros(1, 5), ...
    'cellKm', G.cellKm, 'cacheVersion', 1);
for v = 1:2                                                         % NDSI, NDVI
    fp = fileOf(F, COV{v}, yr);
    checkOrigin(fp, G.origin, G.cell);
    dn = tagDates(fp);  dv = datevec(dn);
    assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fp, yr);
    doy = round(dn - datenum(yr, 1, 0));  S.nDays(v) = numel(doy);
    Fk = openTiff(fp);  ck = onCleanup(@() closeTiff(Fk));
    assert(Fk.H >= G.H1 && Fk.W >= G.W1, '%s is smaller than the Albedo grid', fp);
    for ty = 0:ceil(G.H1 / Fk.th) - 1
        rr = ty*Fk.th + 1 : min((ty+1)*Fk.th, G.H1);
        for tx = 0:ceil(G.W1 / Fk.tw) - 1
            cc = tx*Fk.tw + 1 : min((tx+1)*Fk.tw, G.W1);
            cd = double(reshape(G.code(rr, cc), [], 1));
            if ~any(cd), continue; end
            a = double(readTile(Fk, rr, cc));
            if v == 1, m = a >= 0 & a <= 100; else, m = a >= -2000 & a <= 10000; end
            a(~m) = 0;
            P = sparse(cd + 1, 1:numel(cd), 1, 11, numel(cd));      % region code -> row
            S.sum(:, doy, v) = S.sum(:, doy, v) + full(P * a);
            S.cnt(:, doy, v) = S.cnt(:, doy, v) + full(P * double(m));
        end
    end
    clear ck                                                        % closes the file
end
fe = cellfun(@(x) fileOf(F, x, yr), COV(3:5), 'UniformOutput', false);    % ERA5: t2m, SWE, snowfall
dn = tagDates(fe{1});
for q = 1:3
    checkOrigin(fe{q}, G.origin, G.cell);
    if q > 1, assert(isequal(tagDates(fe{q}), dn), '%s and %s hold different dates', fe{1}, fe{q}); end
end
dv = datevec(dn);
assert(all(dv(:, 1) == yr), '%s: band dates outside %d', fe{1}, yr);
doy = round(dn - datenum(yr, 1, 0));  S.nDays(3:5) = numel(doy);
Ft = openTiff(fe{1});  c1 = onCleanup(@() closeTiff(Ft));
Fw = openTiff(fe{2});  c2 = onCleanup(@() closeTiff(Fw));
Ff = openTiff(fe{3});  c3 = onCleanup(@() closeTiff(Ff));
assert(Fw.th == Ft.th && Ff.th == Ft.th && Fw.tw == Ft.tw && Ff.tw == Ft.tw, 'The ERA5 files are tiled differently');
for ty = 0:ceil(G.H1 / Ft.th) - 1
    rr = ty*Ft.th + 1 : min((ty+1)*Ft.th, G.H1);
    for tx = 0:ceil(G.W1 / Ft.tw) - 1
        cc = tx*Ft.tw + 1 : min((tx+1)*Ft.tw, G.W1);
        cd = double(reshape(G.code(rr, cc), [], 1));
        if ~any(cd), continue; end
        t = double(readTile(Ft, rr, cc));  w = double(readTile(Fw, rr, cc));  f = double(readTile(Ff, rr, cc));
        tok = t ~= -32768 & t ~= 0;
        mw = tok & w ~= -32768 & w ~= 32767;  mf = tok & f ~= -32768;
        t(~tok) = 0;  w(~mw) = 0;  f(~mf) = 0;
        X = {t, w, f};  M = {tok, mw, mf};
        P = sparse(cd + 1, 1:numel(cd), 1, 11, numel(cd));
        for q = 1:3
            S.sum(:, doy, 2 + q) = S.sum(:, doy, 2 + q) + full(P * X{q});
            S.cnt(:, doy, 2 + q) = S.cnt(:, doy, 2 + q) + full(P * double(M{q}));
        end
    end
end
end

% ----------------------------- permafrost zonation -----------------------------

function [pf, pfC] = permafrostGrid(F, G, E)
% Permafrost extent class of every 1-km pixel (pf, H1 x W1) and cell (pfC, HC x WC) from the
% Circum-Arctic Map of Permafrost and Ground-Ice Conditions (Brown et al. 2002, NSIDC GGD318,
% permaice.shp; 1:10,000,000): 4 continuous (90-100 % of the area underlain by permafrost),
% 3 discontinuous (50-90 %), 2 sporadic (10-50 %), 1 isolated patches (0-10 %), 0 land without
% permafrost (and relict permafrost); NaN = glacier, lake, sea, outside the AOI or no polygon.
% A pixel takes the class of the polygon containing its centre. The polygon vertices are moved
% from the map's polar Lambert azimuthal equal-area projection (sphere R = 6,370,997 m, centre the
% North Pole, central meridian 180) to EPSG:3338 (Mapping Toolbox projfwd), then filled per class
% by the even-odd rule over all rings of that class (the map's polygons do not overlap). Cell class
% = the most frequent class of the cell's pixels (ties: the higher class). Cached in E.cacheDir with
% the shapefile's size and date, so a changed file is read again.
fc = fullfile(E.cacheDir, sprintf('permafrost_%dkm.mat', G.cellKm));
d = dir(F.permafrost);
assert(numel(d) == 1, 'Missing the permafrost map %s', F.permafrost);
if isfile(fc)
    L = load(fc);
    if L.bytes == d.bytes && L.fileDate == d.datenum && isequal(size(L.pf), [G.H1 G.W1]) && isequal(size(L.pfC), [G.HC G.WC])
        pf = L.pf;  pfC = L.pfC;  return
    end
end
t0 = tic;
[pd, pn] = fileparts(F.permafrost);
prj = fileread(fullfile(pd, [pn '.prj']));
assert(~isempty(strfind(prj, 'Lambert_Azimuthal_Equal_Area')) && ~isempty(strfind(prj, '6370997')) && ...
    ~isempty(regexp(prj, 'Central_Meridian"?\s*,\s*180', 'once')) && ...
    ~isempty(regexp(prj, 'Latitude_Of_Origin"?\s*,\s*90', 'once')), ...
    ['%s: expected the polar Lambert azimuthal equal-area projection (sphere 6370997 m, central meridian 180) ' ...
     'of the Circum-Arctic permafrost map'], fullfile(pd, [pn '.prj']));
S = shaperead(F.permafrost);
cls = permafrostClass(strtrim({S.EXTENT}), strtrim({S.COMBO}));
p3338 = projcrs(3338);
pf = nan(G.H1, G.W1, 'single');  nOver = 0;
for k = [4 3 2 1 0]
    sel = find(cls == k);
    if isempty(sel), continue; end
    [X, Y] = closeRings([S(sel).X], [S(sel).Y]);
    [lat, lon] = laeaPolarInv(X, Y, 6370997, 180);
    [xa, ya] = projfwd(p3338, lat, lon);
    M = fillRings(xa, ya, G.origin(1), G.origin(2), G.cell, G.H1, G.W1);
    nOver = nOver + nnz(M & isfinite(pf));
    pf(M & isnan(pf)) = k;
end
pf(~G.aoi) = NaN;
pfC = nan(G.HC, G.WC, 'single');  best = zeros(G.HC, G.WC);
for k = 0:4
    n = blockSum(double(pf == k), G);
    pfC(n > 0 & n >= best) = k;  best = max(best, n);
end
bytes = d.bytes;  fileDate = d.datenum;  cellKm = G.cellKm;                        %#ok<NASGU>
save(fc, 'pf', 'pfC', 'bytes', 'fileDate', 'cellKm', '-v7');
fprintf(['Permafrost map: %d polygons, AOI pixels continuous %.1f %%, discontinuous %.1f %%, sporadic %.1f %%, ' ...
    'isolated %.1f %%, none %.1f %%, no class %.1f %% (glacier, lake, sea); %d pixels in two classes (%.0f s)\n'], ...
    numel(S), 100 * mean(pf(G.aoi) == 4), 100 * mean(pf(G.aoi) == 3), 100 * mean(pf(G.aoi) == 2), ...
    100 * mean(pf(G.aoi) == 1), 100 * mean(pf(G.aoi) == 0), 100 * mean(isnan(pf(G.aoi))), nOver, toc(t0));
end

function cls = permafrostClass(extent, combo)
% EXTENT (C, D, S, I) -> 4..1; no extent: COMBO ld (land) and r (relict) -> 0, g (glacier), l (lake)
% and o (sea) -> NaN; anything else is an error
cls = nan(1, numel(extent));
ext = {'I', 'S', 'D', 'C'};
for q = 1:numel(extent)
    k = find(strcmp(ext, extent{q}));
    if ~isempty(k)
        cls(q) = k;
    elseif any(strcmp(combo{q}, {'ld', 'r'}))
        cls(q) = 0;
    else
        assert(any(strcmp(combo{q}, {'g', 'l', 'o'})), 'Unknown permafrost polygon: EXTENT "%s", COMBO "%s"', ...
            extent{q}, combo{q});
    end
end
end

function [X, Y] = closeRings(X, Y)
% NaN-separated rings -> each ring closed (last vertex = first) and followed by NaN
X = X(:)';  Y = Y(:)';
br = [0, find(isnan(X)), numel(X) + 1];
xo = cell(1, numel(br) - 1);  yo = xo;
for q = 1:numel(br) - 1
    i = br(q) + 1 : br(q + 1) - 1;
    if isempty(i), continue; end
    assert(all(isfinite(Y(i))), 'NaN in the Y coordinates of a ring but not in X');
    if X(i(1)) ~= X(i(end)) || Y(i(1)) ~= Y(i(end)), i = [i i(1)]; end
    xo{q} = [X(i) NaN];  yo{q} = [Y(i) NaN];
end
X = [xo{:}];  Y = [yo{:}];
end

function [lat, lon] = laeaPolarInv(x, y, R, lon0)
% inverse of the north-polar Lambert azimuthal equal-area projection on a sphere of radius R
% (Snyder 1987, eqs. 20-14, 20-16): degrees; lon in [-180, 180)
rho = hypot(x, y);
lat = 90 - 2 * asind(min(rho / (2 * R), 1));
lon = mod(lon0 + atan2d(x, -y) + 180, 360) - 180;
end

function M = fillRings(X, Y, x0, y0, cs, H, W)
% even-odd fill of the NaN-separated, closed rings X, Y (m) on the H x W grid with upper-left corner
% (x0, y0) and cell size cs: true where the pixel centre is inside. Pixel (i, j) has its centre at
% (u, v) = (j, i), u = (x - x0) / cs + 0.5, v = (y0 - y) / cs + 0.5. Each row of centres crosses the
% edges with vlo <= i < vhi (half-open, so a vertex on the line counts once); between the 1st and 2nd,
% 3rd and 4th, ... crossing of a row (sorted), the pixels ceil(u1) .. ceil(u2) - 1 are inside.
u = (X(:) - x0) / cs + 0.5;  v = (y0 - Y(:)) / cs + 0.5;
a = [u(1:end - 1), v(1:end - 1)];  b = [u(2:end), v(2:end)];
ok = all(isfinite([a b]), 2) & a(:, 2) ~= b(:, 2);
a = a(ok, :);  b = b(ok, :);
r0 = max(ceil(min(a(:, 2), b(:, 2))), 1);  r1 = min(ceil(max(a(:, 2), b(:, 2))) - 1, H);
n = max(r1 - r0 + 1, 0);  e = repelem((1:numel(n))', n);
if isempty(e), M = false(H, W); return; end
rows = r0(e) + ((1:numel(e))' - repelem(cumsum(n) - n, n)) - 1;
ux = a(e, 1) + (rows - a(e, 2)) .* (b(e, 1) - a(e, 1)) ./ (b(e, 2) - a(e, 2));
Z = sortrows([rows, ux]);
assert(all(mod(accumarray(Z(:, 1), 1, [H 1]), 2) == 0), 'A row crosses the rings an odd number of times (open ring?)');
p1 = Z(1:2:end, :);  p2 = Z(2:2:end, :);
c1 = max(ceil(p1(:, 2)), 1);  c2 = min(ceil(p2(:, 2)) - 1, W);  k = c1 <= c2;
Dm = accumarray([p1(k, 1), c1(k)], 1, [H W + 1]) - accumarray([p1(k, 1), c2(k) + 1], 1, [H W + 1]);
M = cumsum(Dm, 2) > 0;
M = M(:, 1:W);
end

% ----------------------------- stage 2: assemble and analyse -----------------------------

function D = assemble(E, G, smp, COV)
% Cached years and the Prepared layers -> the arrays the analyses use (none of them is
% cells x years: the 8-day and monthly cell values are streamed from the cache later):
%   cnt [H1 x W1]                      valid days over the baseline months (never-observed share)
%   climSum / climN [cells x nM]       sum and number of the baseline years' cell monthly means
%   rsum / rcnt / rmag [11 x 366 x years], nB, valid, mismatch, first, last
%   salb [observed sample pixels x 366 x years], scov.<var> [pixels x days x baseline years]
%   sw.<inAll|nAll> [11 x 366 x baseline years]  incoming shortwave sums (passShortwave)
years = E.years(:)';  nY = numel(years);  ib = years >= E.base(1) & years <= E.base(2);
nM = numel(E.months);  nK = numel(G.keep);  baseYears = years(ib);
D.years = years;  D.ib = ib;  D.months = E.months;
cc = zeros(G.H1, G.W1);
D.climSum = zeros(nK, nM);  D.climN = zeros(nK, nM);
D.rsum = zeros(11, 366, nY);  D.rcnt = D.rsum;  D.rmag = D.rsum;
[D.nB, D.nA, D.nQ, D.valid, D.mismatch, D.first, D.last] = deal(zeros(nY, 1));
salb = zeros(numel(smp.idx), 366, nY, 'int16');
for i = 1:nY
    L = load(fullfile(E.cacheDir, sprintf('alb_%d.mat', years(i))));
    assert(isequal(L.months, E.months) && isequal(size(L.msum), [G.H1 G.W1 nM]) && ...
        size(L.samp, 1) == numel(smp.idx), 'alb_%d.mat was made with other months or another grid: set E.redo = true', years(i));
    if ib(i)
        cm = cellMonthly(L, G, E);                                   % cells x nM
        f = isfinite(cm);  cm(~f) = 0;
        D.climSum = D.climSum + cm;  D.climN = D.climN + f;
        cc = cc + sum(double(L.mcnt), 3);
    end
    D.rsum(:, :, i) = L.rsum;  D.rcnt(:, :, i) = L.rcnt;  D.rmag(:, :, i) = L.rmag;
    D.nB(i) = L.nB;  D.nA(i) = L.nA;  D.nQ(i) = L.nQ;  D.valid(i) = L.valid;  D.mismatch(i) = L.mismatch;
    D.first(i) = L.first;  D.last(i) = L.last;
    salb(:, :, i) = L.samp;
    clear L cm f
end
D.cnt = cc;
% sample pixels with any valid Albedo value
obs = any(reshape(salb ~= -32768, size(salb, 1), []), 2);
D.sIdx = smp.idx(obs);  D.salb = salb(obs, :, :);  D.nSampled = numel(smp.idx);
[r, c] = ind2sub([G.H1 G.W1], D.sIdx);
D.sPixRM = (r - 1) * G.W1 + (c - 1);                             % row-major index (tie-break key)
D.sElev = G.elev(D.sIdx);  D.sSlope = G.slope(D.sIdx);
D.sNorth = northness(G.slope(D.sIdx), G.aspect(D.sIdx));
D.sLat = G.latC(sub2ind([G.HC G.WC], min(floor((r - 1) / G.agg) + 1, G.HC), min(floor((c - 1) / G.agg) + 1, G.WC)));
D.covDays = E.covDays;  dd = E.covDays(1):E.covDays(2);
for v = COV
    D.scov.(v{1}) = zeros(nnz(obs), numel(dd), numel(baseYears), 'int16');
end
fl = {'inAll', 'nAll'};
for f = fl, D.sw.(f{1}) = zeros(11, 366, numel(baseYears)); end
for j = 1:numel(baseYears)
    L = load(fullfile(E.cacheDir, sprintf('sw_%d.mat', baseYears(j))), fl{:});
    for f = fl, D.sw.(f{1})(:, :, j) = L.(f{1}); end
end
for j = 1:numel(baseYears)
    L = load(fullfile(E.cacheDir, sprintf('cov_%d.mat', baseYears(j))));
    for v = COV
        assert(size(L.(v{1}), 1) == numel(smp.idx), 'cov_%d.mat was made for other sample pixels: set E.redo = true', baseYears(j));
        D.scov.(v{1})(:, :, j) = L.(v{1})(obs, dd);
    end
end
end

function S = blockSum(B, G)
% sums over the G.agg x G.agg blocks of the 1-km grid, from the top-left -> HC x WC
a = G.agg;
S = reshape(sum(sum(reshape(B(1:G.HC*a, 1:G.WC*a), a, G.HC, a, G.WC), 1), 3), G.HC, G.WC);
end

function cm = cellMonthly(L, G, E)
% one cached year (msum, mcnt) -> kept cells x nM monthly means: pooled mean of the cell's
% 1-km values, NaN with < E.minCellDays valid pixel-days
nM = numel(E.months);  nK = numel(G.keep);
cs = zeros(nK, nM);  cn = zeros(nK, nM);
ms = double(L.msum);  mc = double(L.mcnt);
for k = 1:nM
    b = blockSum(ms(:, :, k), G);  cs(:, k) = b(G.keep);
    b = blockSum(mc(:, :, k), G);  cn(:, k) = b(G.keep);
end
cm = cs ./ max(cn, 1);  cm(cn < E.minCellDays) = NaN;
end

function B = blockMeans(E, yrs, f, vars)
% 8-day cell means [cells of cache file f x 46 x numel(yrs)] of the variables vars (1..6:
% Albedo, NDSI, NDVI, t2m, SWE, snowfall; raw units, NaN with < E.blkMinCount valid
% pixel-days) -> cell array, one entry per variable
B = cell(1, numel(vars));
for j = 1:numel(yrs)
    L = load(blockFile(E, yrs(j), f), 's', 'n', 'blocks', 'rows', 'cellKm');
    assert(isequal(L.blocks, 1:46) && L.cellKm == E.cellKm, '%s was made for other blocks or cells: set E.redo = true', ...
        blockFile(E, yrs(j), f));
    if j == 1
        nr = L.rows(2) - L.rows(1) + 1;
        for q = 1:numel(vars), B{q} = nan(nr, 46, numel(yrs), 'single'); end
    end
    assert(size(L.s, 1) == nr, '%s holds another number of cells than the other years', blockFile(E, yrs(j), f));
    for q = 1:numel(vars)
        sv = double(L.s(:, :, vars(q)));  nv = double(L.n(:, :, vars(q)));
        mv = sv ./ max(nv, 1);  mv(nv < E.blkMinCount) = NaN;
        B{q}(:, :, j) = single(mv);
    end
end
end

function [val, rclim] = regionalSeries(anom, clim, regC, order, minShare)
% regional value per year (rows: order, then statewide): region mean of the cell
% climatologies + mean anomaly of the cells observed that year (>= minShare of the
% region's cells with a climatology). anom: cells x years; clim: cells x 1.
regs = [order 0];  nY = size(anom, 2);
val = nan(10, nY);  rclim = nan(10, 1);
for j = 1:10
    if regs(j) > 0, inR = regC == regs(j); else, inR = regC > 0; end
    sel = isfinite(clim) & inR;
    if ~any(sel), continue; end
    rclim(j) = mean(clim(sel));
    a = anom(sel, :);
    share = sum(isfinite(a), 1) / nnz(sel);
    v = rclim(j) + nanMean(a, 1);
    v(share < minShare) = NaN;
    val(j, :) = v;
end
end

function R = edaAnalyze(D, G, E, nw)
% nw: parallel workers of the driver maps (one 8-day cache file each); 0 = serial
years = D.years;  ib = D.ib;  nY = numel(years);  nb = nnz(ib);  yb = years(ib);
nM = numel(E.months);  mons = E.months;
code = double(G.code);  npix = G.npix;  nAoi = sum(npix(2:11));
ORDER = G.order;  nK = numel(G.keep);  regC = G.regC;  regs = [ORDER 0];
R.years = years;  R.ib = ib;  R.months = mons;  R.mname = monthNames(mons);
R.base = E.base;  R.seasons = E.seasons;  R.npix = npix;  R.nAoi = nAoi;
R.regNames = G.regNames(ORDER);  R.regShort = G.short(ORDER);       % the 9 regions in plotting order
R.vars = {'Albedo', 'NDSI', 'NDVI', 't2m', 'SWE', 'snowfall', 'elevation', 'slope', 'northness', 'latitude', 'peak NDVI', ...
    'permafrost'};
nV = numel(R.vars);  nStat = 6;                         % the last nStat variables are static (one value per pixel)

% ---------- data by year; never-observed share ----------
aoi = G.aoi;
R.neverShare = nnz(D.cnt == 0 & aoi) / nnz(aoi);
R.regArea = zeros(1, 10);  R.regNever = zeros(1, 10);
for j = 1:10
    if regs(j) > 0, m1 = code == regs(j); else, m1 = code > 0; end
    R.regArea(j) = nnz(m1);  R.regNever(j) = 100 * mean(D.cnt(m1) == 0);
end
stC = reshape(sum(D.rcnt(2:11, :, :), 1), 366, nY);
stM = reshape(sum(D.rmag(2:11, :, :), 1), 366, nY);
yt = struct([]);
for i = 1:nY
    d0 = datenum(years(i), mons(1), 1) - datenum(years(i), 1, 0);
    d1 = datenum(years(i), mons(end) + 1, 1) - 1 - datenum(years(i), 1, 0);
    c = sum(stC(d0:d1, i));  g = sum(stM(d0:d1, i));
    yt(i).year = years(i);  yt(i).bands = D.nB(i);
    yt(i).first = datestr(D.first(i), 'yyyy-mm-dd');  yt(i).last = datestr(D.last(i), 'yyyy-mm-dd');
    yt(i).valid_pixel_days = D.valid(i);
    yt(i).valid_pct = 100 * c / (nAoi * (d1 - d0 + 1));
    yt(i).magnitude_share_pct = 100 * g / c;
    yt(i).qa_mismatch = D.mismatch(i);  yt(i).in_baseline = double(ib(i));
end
R.yearTable = yt;

% ---------- cells: monthly means, climatology, regional series and five-year means (one cached year
% at a time). Period mean of a cell = mean over the months of the period's mean monthly value (a month
% needs >= E.perMinFrac of the period's years; every month needed), so that every period weighs the
% months alike; R.perBase = the same from the baseline climatology ----------
cclim = D.climSum ./ D.climN;  cclim(D.climN < E.minClimYears) = NaN;    % cells x nM
R.regVal = nan(10, nY, nM);  R.regClim = nan(10, nM);
Pp = E.periods;  nP = size(Pp, 1);
assert(size(Pp, 2) == 2 && all(Pp(:, 1) <= Pp(:, 2)) && all(Pp(2:end, 1) > Pp(1:end - 1, 2)) && issorted(years), ...
    'E.periods must be [first last] rows in increasing order without overlap');
R.periods = Pp;  R.perYears = zeros(nP, 1);  R.perMean = nan(nK, nP, 'single');
pS = zeros(nK, nM);  pN = zeros(nK, nM);
for i = 1:nY
    L = load(fullfile(E.cacheDir, sprintf('alb_%d.mat', years(i))), 'msum', 'mcnt');
    cm = cellMonthly(L, G, E);                                           % cells x nM
    for k = 1:nM
        [R.regVal(:, i, k), R.regClim(:, k)] = regionalSeries(cm(:, k) - cclim(:, k), cclim(:, k), regC, ORDER, E.minRegShare);
    end
    q = find(years(i) >= Pp(:, 1) & years(i) <= Pp(:, 2));
    if isempty(q), continue; end
    f = isfinite(cm);  cm(~f) = 0;  pS = pS + cm;  pN = pN + f;
    if years(i) == Pp(q, 2) || i == nY || years(i + 1) > Pp(q, 2)         % last year of the period read
        ny = nnz(years >= Pp(q, 1) & years <= Pp(q, 2));
        pm = pS ./ max(pN, 1);  pm(pN < ceil(E.perMinFrac * ny - 1e-9)) = NaN;
        R.perMean(:, q) = single(mean(pm, 2));                               % NaN unless every month
        R.perYears(q) = ny;  pS(:) = 0;  pN(:) = 0;
    end
end
assert(all(R.perYears > 0), 'A period in E.periods holds none of the years read (E.years)');
R.perBase = single(mean(cclim, 2));
common = all(isfinite(R.perMean), 2) & isfinite(R.perBase);              % cells with every period
R.perReg = nan(10, nP);  R.perRegBase = nan(10, 1);  R.perRegN = zeros(10, 1);
for j = 1:10
    if regs(j) > 0, sel = common & regC == regs(j); else, sel = common; end
    R.perRegN(j) = nnz(sel);
    if any(sel), R.perReg(j, :) = mean(double(R.perMean(sel, :)), 1);  R.perRegBase(j) = mean(double(R.perBase(sel))); end
end
clear L cm cclim pS pN pm f common

% ---------- mean annual cycle (daily, regions and statewide) ----------
npx = [npix(ORDER + 1), nAoi];
daily = nan(10, 365, nY);
for i = 1:nY
    Sd = to365([D.rsum(ORDER + 1, :, i); sum(D.rsum(2:11, :, i), 1)], years(i));
    Cd = to365([D.rcnt(ORDER + 1, :, i); sum(D.rcnt(2:11, :, i), 1)], years(i));
    v = Sd ./ max(Cd, 1);  v(Cd < E.minRegDayCov * npx(:)) = NaN;
    daily(:, :, i) = v;
end
R.cycMean = nan(10, 365);
for j = 1:10
    X = reshape(daily(j, :, ib), 365, nb);
    n = sum(isfinite(X), 2);
    mu = nanMean(X, 2);  mu(n < E.minCycleYears) = NaN;
    R.cycMean(j, :) = nanMovMean(mu, 3, 4);
end
ibIdx = find(ib);  R.cycYear = nan(10, 365, nb, 'single');             % every baseline year, 7-day smoothed
for i = 1:nb
    for j = 1:10, R.cycYear(j, :, i) = single(nanMovMean(daily(j, :, ibIdx(i)), 3, 4)); end
end

% ---------- incoming shortwave: daily climatology of the ERA5-Land SSRD (regions and statewide; stage 2h) ----------
R.swIn = shortwave(D.sw, yb, E, ORDER);

% ---------- 8-day blocks: driver-anomaly correlations per cell ----------
% anomaly = 8-day cell mean - the cell's mean of that block over the baseline; per cell and
% season: Spearman rho and Chatterjee xi(driver -> Albedo) over the (year, block) pairs.
% One cache file of cells at a time (all baseline years), E.blkChunk cells per step.
drv = {'t2m', 'swe', 'snowfall', 'ndsi', 'ndvi'};  dv = [4 5 6 2 3];   % their variable in the 8-day cache
nS = size(E.blkSeasons, 1);  nD = numel(drv);  blocks = 1:46;
[rC, cC] = ind2sub([G.HC G.WC], G.keep);
cellRM = (rC(:) - 1) * G.WC + (cC(:) - 1);                            % row-major cell index (tie key)
R.blkRho = nan(nS, nD, nK);  R.blkXi = nan(nS, nD, nK);  R.blkN = zeros(nS, nD, nK);
R.blkSeasons = E.blkSeasons(:, 1)';  R.blkDrivers = drv;
nF = ceil(nK / E.blkFileCells);  outF = cell(1, nF);  t1 = tic;
try
    parfor (f = 1:nF, nw)
        outF{f} = driverMapsFile(E, yb, f, nK, cellRM, dv);
    end
catch err
    if nw == 0, rethrow(err); end
    fprintf('  the parallel driver maps failed (%s): running them serially\n', err.message);
    for f = 1:nF, outF{f} = driverMapsFile(E, yb, f, nK, cellRM, dv); end
end
for f = 1:nF
    rows = outF{f}.rows;
    R.blkRho(:, :, rows) = outF{f}.rho;  R.blkXi(:, :, rows) = outF{f}.xi;  R.blkN(:, :, rows) = outF{f}.n;
    outF{f} = [];
end
fprintf('  driver maps: %d cells (%.0f s)\n', nK, toc(t1));
clear outF
% the driver-anomaly correlations grouped by the permafrost extent class of the cell (static: it has
% no anomaly of its own): per season, driver and class the quartiles of rho and the median xi
R.pfZones = {'none', 'isolated', 'sporadic', 'discontinuous', 'continuous'};   % classes 0..4
pfK = double(G.pfC(G.keep));
R.pfZoneRho = nan(nS, nD, 5, 3);  R.pfZoneXi = nan(nS, nD, 5);  R.pfZoneN = zeros(nS, nD, 5);
for si = 1:nS
    for di = 1:nD
        r = squeeze(R.blkRho(si, di, :));  x = squeeze(R.blkXi(si, di, :));
        for k = 0:4
            sel = pfK == k & isfinite(r);  R.pfZoneN(si, di, k + 1) = nnz(sel);
            if any(sel)
                R.pfZoneRho(si, di, k + 1, :) = [pctl(r(sel), 25), median(r(sel)), pctl(r(sel), 75)];
                R.pfZoneXi(si, di, k + 1) = median(x(sel));
            end
        end
    end
end
R.pfCells = arrayfun(@(k) nnz(pfK == k), 0:4);

% ---------- drivers: Spearman rho and Chatterjee xi (sample pixels, baseline) ----------
nd = size(D.scov.ndsi, 2);  days = D.covDays(1) + (0:nd - 1);           % actual day of year
albRaw = D.salb(:, days, ib);                                           % pixels x days x years
nvd = double(D.scov.ndvi(:, days >= 182 & days <= 243, :));            % Jul-Aug
nvd(nvd < -2000 | nvd > 10000) = NaN;
pkY = reshape(max(nvd, [], 2), size(nvd, 1), []) * 1e-4;               % pixels x years
peak = nanMean(pkY, 2);  peak(sum(isfinite(pkY), 2) < E.peakMinYears) = NaN;
clear nvd pkY
stat = [D.sElev(:), D.sSlope(:), D.sNorth(:), D.sLat(:), peak(:), double(G.pf(D.sIdx(:)))];   % permafrost: class 0..4
nP = numel(D.sIdx);  pkey = tieKey(D.sPixRM(:), 0, 0);
R.corr = cell(1, 3);  R.corrN = cell(1, 3);  R.rho = cell(1, 3);  R.xi = cell(1, 3);  R.xiN = cell(1, 3);
for s = 1:3
    dsel = find(days >= E.seasons{s, 2} & days <= E.seasons{s, 3});
    sub = albRaw(:, dsel, :);
    lin = find(sub ~= -32768);                                          % pixel-days with an Albedo value
    [pa, dd, ya] = ind2sub(size(sub), lin);
    dcol = reshape(dsel(dd), [], 1);
    li = sub2ind(size(D.scov.ndsi), pa, dcol, ya);
    V = nan(numel(lin), nV);
    V(:, 1) = double(sub(lin)) * 1e-4;
    x = double(D.scov.ndsi(li));  x(x < 0 | x > 100) = NaN;  V(:, 2) = x;
    x = double(D.scov.ndvi(li));  x(x < -2000 | x > 10000) = NaN;  V(:, 3) = x * 1e-4;
    x = double(D.scov.t2m(li));  tok = x ~= -32768 & x ~= 0;
    V(:, 4) = x * 0.01 - 273.15;  V(~tok, 4) = NaN;
    x = double(D.scov.swe(li));  V(:, 5) = x * 0.1;  V(~(tok & x ~= -32768 & x ~= 32767), 5) = NaN;
    x = double(D.scov.snowfall(li));  V(:, 6) = x * 0.1;  V(~(tok & x ~= -32768), 6) = NaN;
    V(:, nV - nStat + 1:nV) = stat(pa, :);
    key = tieKey(D.sPixRM(pa), yb(ya)', days(dcol)');
    clear sub lin dd li x tok
    [R.corr{s}, R.corrN{s}] = spearmanMatrix(V);
    rho = nan(1, nV);  xi = nan(1, nV);  xn = zeros(1, nV);
    for j = 2:6                                                         % daily covariates: pixel-days
        m = isfinite(V(:, j)) & isfinite(V(:, 1));  xn(j) = nnz(m);
        if xn(j) > 2, rho(j) = R.corr{s}(1, j);  xi(j) = xiCorr(V(m, j), V(m, 1), key(m)); end
    end
    cnt = accumarray(pa, 1, [nP 1]);                                    % static covariates: pixels
    pm = accumarray(pa, V(:, 1), [nP 1]) ./ max(cnt, 1);  pm(cnt < E.minPixDays) = NaN;
    for j = nV - nStat + 1:nV
        st = stat(:, j - nV + nStat);  m = isfinite(pm) & isfinite(st);  xn(j) = nnz(m);
        if xn(j) > 2
            rho(j) = pearson(rankavg(st(m)), rankavg(pm(m)));  xi(j) = xiCorr(st(m), pm(m), pkey(m));
        end
    end
    R.rho{s} = rho;  R.xi{s} = xi;  R.xiN{s} = xn;
    clear V pa ya key
end
R.nSampled = D.nSampled;  R.nObserved = numel(D.sIdx);  R.sIdx = D.sIdx;
end

function o = driverMapsFile(E, yb, f, nK, cellRM, dv)
% driver maps of the cells of 8-day cache file f (one parfor iteration): per season and driver, Spearman rho,
% Chatterjee xi and the pairs [nS x nD x rows] (E.blkChunk cells at a time, which does not change the results)
nS = size(E.blkSeasons, 1);  nD = numel(dv);  blocks = 1:46;  nb = numel(yb);
rows = (f - 1) * E.blkFileCells + 1 : min(f * E.blkFileCells, nK);  nr = numel(rows);
B = blockMeans(E, yb, f, 1:6);                                    % Albedo and drivers of these cells
o.rows = rows;
o.rho = nan(nS, nD, nr);  o.xi = nan(nS, nD, nr);  o.n = zeros(nS, nD, nr);
for c0 = 1:E.blkChunk:nr
    loc = c0:min(c0 + E.blkChunk - 1, nr);  cs = rows(loc);  nc_ = numel(cs);
    for si = 1:nS
        cols = find(ismember(blocks, E.blkSeasons{si, 2}));  nc = numel(cols);
        key = tieKey(repmat(cellRM(cs), 1, nc * nb), repmat(kron(yb, ones(1, nc)), nc_, 1), ...
                     repmat(repmat(blocks(cols), 1, nb), nc_, 1));
        A = blockAnomalies(B{1}(loc, cols, :), E.blkMinYears);
        for di = 1:nD
            X = blockAnomalies(B{dv(di)}(loc, cols, :), E.blkMinYears);
            [r, n] = rowSpearman(X, A, E.blkMinPairs);
            o.rho(si, di, loc) = r;  o.n(si, di, loc) = n;
            o.xi(si, di, loc) = rowXi(X, A, key, E.blkMinPairs);
        end
    end
end
fprintf('  driver maps: cells %d-%d of %d\n', rows(1), rows(end), nK);
end

% ----------------------------- driver cycles, annual metrics, dependence by region, snow years -----------------------------

function X = edaExtra(R, G, E, smp, COV)
% Stage 2b, from the cache and R (rows of every [10 x ...] array: the 9 ecoregions in G.order, then statewide, as
% R.cycMean):
%   drvCycYear [10 x 365 x baseline years x 5], drvCycMean [10 x 365 x 5]  annual cycle of each driver, as the
%       Albedo cycle of fig01 (pooled regional daily mean over the valid pixels, >= E.minRegDayCov of them; mean
%       over >= E.minCycleYears years; 7-day centred mean), in physical units
%   met [10 x baseline years x nMet], metTrend  annual Albedo metrics of every region and year from the year's
%       7-day smoothed regional curve (R.cycYear) and their Theil-Sen trends (cycleMetrics, senKendall)
%   dep  Spearman rho and Chatterjee xi of Albedo vs each covariate by region, as fig03 (dependenceByRegion)
%   snow statewide daily Albedo and NDSI snow cover of the E.snowYears windows (snowYearSeries)
yb = R.years(R.ib);  nb = numel(yb);  ORDER = G.order;
npx = [G.npix(ORDER + 1), sum(G.npix(2:11))];
X.yb = yb;  X.base = E.base;
% ---------- driver annual cycles ----------
X.drv = COV;
X.drvName  = {'NDSI snow cover', 'NDVI', '2-m air temperature', 'snow water equivalent', 'snowfall'};
X.drvUnit  = {'0-1', '', 'deg C', 'mm', 'mm d-1'};
X.drvScale = [0.01 1e-4 0.01 0.1 0.1];  X.drvOffset = [0 0 -273.15 0 0];    % raw -> physical
daily = nan(10, 365, nb, 5);
for i = 1:nb
    L = loadDrv(E, yb(i), COV);
    for v = 1:5
        Sd = to365([L.sum(ORDER + 1, :, v); sum(L.sum(2:11, :, v), 1)], yb(i));
        Cd = to365([L.cnt(ORDER + 1, :, v); sum(L.cnt(2:11, :, v), 1)], yb(i));
        x = Sd ./ max(Cd, 1) * X.drvScale(v) + X.drvOffset(v);  x(Cd < E.minRegDayCov * npx(:)) = NaN;
        daily(:, :, i, v) = x;
    end
end
X.drvCycMean = nan(10, 365, 5);  X.drvCycYear = nan(10, 365, nb, 5, 'single');
for v = 1:5
    for j = 1:10
        Y = reshape(daily(j, :, :, v), 365, nb);
        n = sum(isfinite(Y), 2);  mu = nanMean(Y, 2);  mu(n < E.minCycleYears) = NaN;
        X.drvCycMean(j, :, v) = nanMovMean(mu, 3, 4);
        for i = 1:nb, X.drvCycYear(j, :, i, v) = single(nanMovMean(Y(:, i), 3, 4)); end
    end
end
clear daily
% ---------- annual Albedo metrics of every region and year, and their trends ----------
[X.met, X.metName, X.metUnit] = deal([]);
for j = 1:10
    for i = 1:nb
        [m, X.metName, X.metUnit] = cycleMetrics(double(R.cycYear(j, :, i)), R.cycMean(j, :), E.metricMinFrac);
        if isempty(X.met), X.met = nan(10, nb, numel(m)); end
        X.met(j, i, :) = m;
    end
end
nMet = size(X.met, 3);
X.metTrend = struct('slope', nan(10, nMet), 'tau', nan(10, nMet), 'p', nan(10, nMet), 'n', zeros(10, nMet));
for j = 1:10
    for q = 1:nMet
        [X.metTrend.slope(j, q), X.metTrend.tau(j, q), X.metTrend.p(j, q), X.metTrend.n(j, q)] = ...
            senKendall(yb, X.met(j, :, q), 10);
    end
end
X.metTrend.q = bhFdr(X.metTrend.p);             % Benjamini-Hochberg over all region x metric tests
% ---------- Spearman rho and Chatterjee xi of Albedo vs the covariates, by region ----------
X.dep = dependenceByRegion(G, E, smp, COV, yb);
X.dep.checkStatewide = 0;                       % statewide column vs R.rho / R.xi of the main analysis
for s = 1:3
    d = [X.dep.rho(10, :, s) - R.rho{s}(2:end), X.dep.xi(10, :, s) - R.xi{s}(2:end)];
    X.dep.checkStatewide = max([X.dep.checkStatewide, abs(d(isfinite(d)))]);
end
% ---------- statewide Albedo and NDSI snow cover of the snow-year windows ----------
X.snow = snowYearSeries(R, G, E, X, COV);
end

function S = extraSettings(E, R)
% what edaExtra's results depend on: the settings it reads, the sample pixels' seed and size, and a stamp of
% the main results R (eda_results.mat). eda_extra.mat is reused (E.figuresOnly) only when these are equal;
% version: raise it whenever edaExtra or its local functions change
S = struct('version', 2, 'cellKm', E.cellKm, 'years', E.years, 'base', E.base, 'seasons', {E.seasons}, ...
    'covDays', E.covDays, 'minRegDayCov', E.minRegDayCov, 'minCycleYears', E.minCycleYears, ...
    'metricMinFrac', E.metricMinFrac, 'minPixDays', E.minPixDays, 'peakMinYears', E.peakMinYears, ...
    'regMinN', E.regMinN, 'snowYears', E.snowYears, 'nSample', E.nSample, 'seed', E.seed, ...
    'rStamp', [nnz(isfinite(R.cycMean)), sum(R.cycMean(isfinite(R.cycMean))), sum(cellfun(@(r) sum(r(isfinite(r))), R.rho))]);
end

function L = loadDrv(E, yr, COV)
L = load(drvFile(E, yr));
assert(L.cellKm == E.cellKm && isequal(L.vars, COV) && isequal(size(L.sum), [11 366 5]), ...
    '%s was made for other cells or variables: set E.redo = true', drvFile(E, yr));
end

function [m, name, unit] = cycleMetrics(c, clim, minFrac)
% Annual metrics of one year's regional Albedo curve c (1 x 365, 7-day smoothed, 365-day calendar; NaN =
% no value); clim (1 x 365) is the region's baseline mean curve (R.cycMean). A window needs >= minFrac of its
% days with a value, otherwise its metrics are NaN.
%   1 mean Mar-Oct (days 60-304)       2 spring peak: maximum over Mar-Jun (60-181)   3 its day of year
%   4 summer mean Jun-Aug (152-243)    5 summer minimum (152-243)                   6 amplitude = peak - minimum
%   7 half-melt date: first day after the peak below peak - amplitude / 2 (linear between the days with a value)
%   8 melt duration: days from peak - 0.2 amplitude to peak - 0.8 amplitude (after the peak)
%   9 autumn half-return date: first day after the summer minimum (to 31 Oct) at or above minimum + amplitude / 2
%  10 snow-free season = autumn half-return - half-melt (days)
% The window means (1, 4) are the climatological mean of the window plus the year's mean anomaly from clim
% over its days with a value: missing days (low sun in March, cloud) then do not shift the mean towards the
% part of the season that was observed, so years with different gaps stay comparable (window = the days
% with a climatology). The peak date (3) is NaN when the maximum falls on the window's first day with a value:
% the curve may have been higher before (melt already under way on 1 Mar, or no retrieval earlier), so the
% date is censored; the peak value is then a lower bound and is kept.
name = {'mean Mar-Oct', 'spring peak', 'spring peak date', 'summer mean JJA', 'summer minimum', 'amplitude', ...
        'half-melt date', 'melt duration 20-80 %', 'autumn half-return date', 'snow-free season'};
unit = {'Albedo', 'Albedo', 'day of year', 'Albedo', 'Albedo', 'Albedo', 'day of year', 'days', 'day of year', 'days'};
assert(isequal(size(c), [1 365]) && isequal(size(clim), [1 365]), 'cycleMetrics expects 1 x 365 curves');
m = nan(1, 10);
ok = @(d) mean(isfinite(c(d))) >= minFrac;
w = 60:304;  w = w(isfinite(clim(w)));
if ~isempty(w) && ok(w), v = w(isfinite(c(w)));  m(1) = mean(clim(w)) + mean(c(v) - clim(v)); end
ws = 60:181;  wm = 152:243;  wmc = wm(isfinite(clim(wm)));
dPk = NaN;
if ok(ws)
    [m(2), k] = max(c(ws));  dPk = ws(k);                            % max ignores NaN
    if dPk > ws(find(isfinite(c(ws)), 1)), m(3) = dPk; end           % else censored at the window start
end
if ok(wm)
    [m(5), k] = min(c(wm));  dMin = wm(k);
    if ~isempty(wmc) && ok(wmc), v = wmc(isfinite(c(wmc)));  m(4) = mean(clim(wmc)) + mean(c(v) - clim(v)); end
end
if ~(isfinite(m(2)) && isfinite(m(5))), return; end
m(6) = m(2) - m(5);
m(7) = crossDay(c, dPk, 304, m(2) - 0.5 * m(6), -1);
m(8) = crossDay(c, dPk, 304, m(2) - 0.8 * m(6), -1) - crossDay(c, dPk, 304, m(2) - 0.2 * m(6), -1);
m(9) = crossDay(c, dMin, 304, m(5) + 0.5 * m(6), 1);
m(10) = m(9) - m(7);
end

function d = crossDay(c, d0, d1, h, dirn)
% first day after d0 (to d1) at which c falls below h (dirn = -1) or rises to >= h (dirn = +1), linear
% between the neighbouring days with a value (a gap is bridged); NaN if it does not happen
dd = d0:d1;  dd = dd(isfinite(c(dd)));  d = NaN;
for k = 2:numel(dd)
    if (dirn < 0 && c(dd(k)) < h) || (dirn > 0 && c(dd(k)) >= h)
        a = dd(k - 1);  b = dd(k);
        d = a + (b - a) * (c(a) - h) / (c(a) - c(b));
        return
    end
end
end

function [slope, tau, p, n, S, v] = senKendall(t, y, minN)
% Theil-Sen slope per decade and Mann-Kendall tau with its two-sided p (normal approximation, variance
% corrected for ties; no correction for autocorrelation) over the finite values; NaN with < minN values.
% S: the Mann-Kendall score, v: its variance (for mkPCorrected, stage 2g)
ok = isfinite(y(:));  t = t(:);  t = t(ok);  y = y(:);  y = y(ok);  n = numel(y);
slope = NaN;  tau = NaN;  p = NaN;  S = NaN;  v = NaN;
if n < minN, return; end
[i, j] = find(triu(true(n), 1));
dy = y(j) - y(i);  dt = t(j) - t(i);
slope = 10 * median(dy ./ dt);
S = sum(sign(dy) .* sign(dt));
[~, ~, g] = unique(y);  tp = accumarray(g, 1);
v = (n * (n - 1) * (2 * n + 5) - sum(tp .* (tp - 1) .* (2 * tp + 5))) / 18;
tau = S / (n * (n - 1) / 2);
z = (S - sign(S)) / sqrt(v);  if S == 0, z = 0; end
p = erfc(abs(z) / sqrt(2));
end

function q = bhFdr(p)
% Benjamini-Hochberg adjusted p (q) of the finite p in the array p (same size; NaN stays NaN)
q = nan(size(p));  k = find(isfinite(p));  n = numel(k);
if n == 0, return; end
[ps, o] = sort(p(k));
qs = min(1, ps(:) .* n ./ (1:n)');
qs = flipud(cummin(flipud(qs)));
q(k(o)) = qs;
end

function dep = dependenceByRegion(G, E, smp, COV, yb)
% Spearman rho and Chatterjee xi(covariate -> Albedo) of fig03 computed within each ecoregion (and
% statewide): the sample pixels of the region (region code of the pixel's cell), with exactly the variables,
% pixel-days, static covariates and tie-break keys of the statewide analysis in edaAnalyze. A value needs
% >= E.regMinN pixel-days (daily covariates) or pixels (static ones). dep.rho / xi / n [10 x 11 x 3]:
% regions (G.order, then statewide) x covariates (R.vars(2:end)) x seasons.
years = E.years;  nY = numel(years);  ib = years >= E.base(1) & years <= E.base(2);
assert(isequal(years(ib), yb), 'The baseline years differ from those of the results');
salb = zeros(numel(smp.idx), 366, nY, 'int16');
for i = 1:nY
    L = load(fullfile(E.cacheDir, sprintf('alb_%d.mat', years(i))), 'samp');
    assert(size(L.samp, 1) == numel(smp.idx), 'alb_%d.mat was made for other sample pixels', years(i));
    salb(:, :, i) = L.samp;
end
obs = any(reshape(salb ~= -32768, size(salb, 1), []), 2);        % as assemble: observed in any year read
sIdx = smp.idx(obs);  salb = salb(obs, :, ib);
[r, c] = ind2sub([G.H1 G.W1], sIdx);
sPixRM = (r - 1) * G.W1 + (c - 1);
sLat = G.latC(sub2ind([G.HC G.WC], min(floor((r - 1) / G.agg) + 1, G.HC), min(floor((c - 1) / G.agg) + 1, G.WC)));
dd = E.covDays(1):E.covDays(2);
for v = COV, scov.(v{1}) = zeros(nnz(obs), numel(dd), numel(yb), 'int16'); end
for j = 1:numel(yb)
    L = load(fullfile(E.cacheDir, sprintf('cov_%d.mat', yb(j))));
    for v = COV, scov.(v{1})(:, :, j) = L.(v{1})(obs, dd); end
end
days = dd;  nV = 12;  nStat = 6;
albRaw = salb(:, days, :);
nvd = double(scov.ndvi(:, days >= 182 & days <= 243, :));  nvd(nvd < -2000 | nvd > 10000) = NaN;
pkY = reshape(max(nvd, [], 2), size(nvd, 1), []) * 1e-4;
peak = nanMean(pkY, 2);  peak(sum(isfinite(pkY), 2) < E.peakMinYears) = NaN;
stat = [G.elev(sIdx), G.slope(sIdx), northness(G.slope(sIdx), G.aspect(sIdx)), sLat(:), peak(:), double(G.pf(sIdx))];
pkey = tieKey(sPixRM, 0, 0);  nP = numel(sIdx);
pcode = double(G.code(sIdx));  pcode = pcode(:);
inReg = false(nP, 10);
for j = 1:9, inReg(:, j) = pcode == G.order(j); end
inReg(:, 10) = pcode > 0;
dep.rho = nan(10, nV - 1, 3);  dep.xi = dep.rho;  dep.n = zeros(10, nV - 1, 3);  dep.npix = sum(inReg, 1);
for s = 1:3
    dsel = find(days >= E.seasons{s, 2} & days <= E.seasons{s, 3});
    sub = albRaw(:, dsel, :);
    lin = find(sub ~= -32768);
    [pa, di, ya] = ind2sub(size(sub), lin);
    dcol = reshape(dsel(di), [], 1);
    li = sub2ind(size(scov.ndsi), pa, dcol, ya);
    V = nan(numel(lin), nV);
    V(:, 1) = double(sub(lin)) * 1e-4;
    x = double(scov.ndsi(li));  x(x < 0 | x > 100) = NaN;  V(:, 2) = x;
    x = double(scov.ndvi(li));  x(x < -2000 | x > 10000) = NaN;  V(:, 3) = x * 1e-4;
    x = double(scov.t2m(li));  tok = x ~= -32768 & x ~= 0;
    V(:, 4) = x * 0.01 - 273.15;  V(~tok, 4) = NaN;
    x = double(scov.swe(li));  V(:, 5) = x * 0.1;  V(~(tok & x ~= -32768 & x ~= 32767), 5) = NaN;
    x = double(scov.snowfall(li));  V(:, 6) = x * 0.1;  V(~(tok & x ~= -32768), 6) = NaN;
    V(:, nV - nStat + 1:nV) = stat(pa, :);
    key = tieKey(sPixRM(pa), yb(ya)', days(dcol)');
    clear sub lin di li x tok
    cnt = accumarray(pa, 1, [nP 1]);
    pm = accumarray(pa, V(:, 1), [nP 1]) ./ max(cnt, 1);  pm(cnt < E.minPixDays) = NaN;
    for j = 1:10
        rowj = inReg(pa, j);
        for q = 2:6                                                  % daily covariates: pixel-days
            m = rowj & isfinite(V(:, q)) & isfinite(V(:, 1));  dep.n(j, q - 1, s) = nnz(m);
            if nnz(m) >= max(E.regMinN, 3)
                dep.rho(j, q - 1, s) = pearson(rankavg(V(m, q)), rankavg(V(m, 1)));
                dep.xi(j, q - 1, s) = xiCorr(V(m, q), V(m, 1), key(m));
            end
        end
        for q = nV - nStat + 1:nV                                    % static covariates: pixels
            st = stat(:, q - nV + nStat);  m = inReg(:, j) & isfinite(pm) & isfinite(st);  dep.n(j, q - 1, s) = nnz(m);
            if nnz(m) >= max(E.regMinN, 3)
                dep.rho(j, q - 1, s) = pearson(rankavg(st(m)), rankavg(pm(m)));
                dep.xi(j, q - 1, s) = xiCorr(st(m), pm(m), pkey(m));
            end
        end
    end
    clear V pa ya key
end
end

function S = snowYearSeries(R, G, E, X, COV)
% statewide pooled daily mean Albedo and NDSI snow cover (0-1) from 1 Sep of each year in E.snowYears to
% 30 Aug of the next (a day needs >= E.minRegDayCov of the AOI pixels valid), their 7-day centred means and
% the baseline climatologies of the same calendar days (R.cycMean, X.drvCycMean; none on 29 Feb)
S = struct('y0', {}, 'dn', {}, 'alb', {}, 'ndsi', {}, 'albS', {}, 'ndsiS', {}, 'albClim', {}, 'ndsiClim', {}, ...
    'albCover', {}, 'ndsiCover', {});
nAoi = sum(G.npix(2:11));
for k = 1:numel(E.snowYears)
    y0 = E.snowYears(k);  dn = (datenum(y0, 9, 1):datenum(y0 + 1, 8, 30))';
    dv = datevec(dn);  doy = dn - datenum(dv(:, 1), 1, 0);
    leap = (mod(dv(:, 1), 4) == 0 & mod(dv(:, 1), 100) ~= 0) | mod(dv(:, 1), 400) == 0;
    d365 = doy - (leap & doy > 60);  d365(leap & doy == 60) = NaN;
    a = nan(numel(dn), 1);  s = a;  ca = a;  cs = a;
    for y = [y0, y0 + 1]
        sel = dv(:, 1) == y;
        A = load(fullfile(E.cacheDir, sprintf('alb_%d.mat', y)), 'rsum', 'rcnt');
        Dd = loadDrv(E, y, COV);
        sa = sum(A.rsum(2:11, doy(sel)), 1)';  na = sum(A.rcnt(2:11, doy(sel)), 1)';
        sn = sum(Dd.sum(2:11, doy(sel), 1), 1)';  nn = sum(Dd.cnt(2:11, doy(sel), 1), 1)';
        x = sa ./ max(na, 1);  x(na < E.minRegDayCov * nAoi) = NaN;  a(sel) = x;  ca(sel) = na / nAoi;
        x = sn ./ max(nn, 1) * X.drvScale(1);  x(nn < E.minRegDayCov * nAoi) = NaN;  s(sel) = x;  cs(sel) = nn / nAoi;
    end
    ac = nan(numel(dn), 1);  sc = ac;  ok = isfinite(d365);
    ac(ok) = R.cycMean(10, d365(ok));  sc(ok) = X.drvCycMean(10, d365(ok), 1);
    S(k) = struct('y0', y0, 'dn', dn, 'alb', a, 'ndsi', s, 'albS', nanMovMean(a, 3, 4)', 'ndsiS', nanMovMean(s, 3, 4)', ...
        'albClim', ac, 'ndsiClim', sc, 'albCover', ca, 'ndsiCover', cs);
end
end

function printExtra(X, E)
nb = numel(X.yb);  q = [1 2 7 10];                                     % metrics printed
for k = q
    fprintf('Albedo %s, statewide: Theil-Sen %+.3g %s per decade, Kendall tau %+.2f (p %.2f, FDR q %.2f, %d years)\n', ...
        X.metName{k}, X.metTrend.slope(10, k), X.metUnit{k}, X.metTrend.tau(10, k), X.metTrend.p(10, k), ...
        X.metTrend.q(10, k), X.metTrend.n(10, k));
end
nOk = squeeze(sum(isfinite(X.met), 2));                               % [10 x nMet] years with a value
fprintf('Years with a value per metric (min over regions): %s\n', strjoin(arrayfun(@(k) sprintf('%s %d', ...
    X.metName{k}, min(nOk(:, k))), 1:size(nOk, 2), 'UniformOutput', false), ', '));
fprintf('Annual metrics: %d regions x %d years; fig03 by region: max |statewide - main analysis| = %.3g\n', 10, nb, ...
    X.dep.checkStatewide);
for k = 1:numel(X.snow)
    sw = X.snow(k);
    fprintf('Snow year %d-%d: %d days with statewide Albedo, %d with NDSI\n', sw.y0, sw.y0 + 1, nnz(isfinite(sw.alb)), ...
        nnz(isfinite(sw.ndsi)));
end
if isempty(E.snowYears), fprintf('No snow-year figure (E.snowYears is empty)\n'); end
end

% ----------------------------- durations: Albedo, NDSI snow cover and air temperature (stage 2g) ----------------------

function DU = edaDurations(R, G, E, X)
% Stage 2g: how many days each year Albedo, NDSI snow cover and 2-m air temperature stay past a threshold (a season) and
% near their summer extreme (at peak), in every region (rows as R.cycMean: the 9 ecoregions in G.order, then All
% Alaska), from each baseline year's 7-day smoothed regional curves (those of fig01: Albedo R.cycYear, NDSI and t2m
% X.drvCycYear):
%   1 Albedo season       Albedo below E.durAlbedo
%   2 Albedo at peak      Albedo below E.durAlbedo and less than E.durPeakAlbedo above the extreme
%   3 snow-free season    NDSI snow cover below E.durNdsi
%   4 NDSI at peak        NDSI below E.durNdsi and less than E.durPeakNdsi above the extreme
%   5 thaw season         t2m above E.durT2m
%   6 warmest weeks       t2m above E.durT2m and less than E.durPeakT2m below the extreme
% The extreme (the peak) is the year's most summer-like value of the curve in the days E.durPivot (lowest Albedo and
% NDSI, highest t2m), its day the pivot; the winter extremes cannot be measured (no Albedo in the polar night, no NDSI
% from October to February). Gaps of <= E.durMaxGap days are filled linearly (fillGaps). A season is the run of days
% around the pivot from the day after the last spell of >= E.durPersist consecutive days out of its state before the
% pivot to the day before the first such spell after it (seasonRun): shorter spells inside belong to it. A search that
% meets a longer gap or the end of the year leaves that boundary unknown (NaN); a pivot out of the state (the season's
% threshold not passed) gives 0 days, and then 0 days at peak. The days at peak are every day of the year in the peak's
% state (stateDays; a run around the pivot would break at a single noisy week), unknown (NaN) when one of them touches a
% day without a value or the start or end of the year. A curve with < E.metricMinFrac of the days of E.durPivot observed
% gives no values that year. Trends of the six durations: Theil-Sen per decade, Mann-Kendall (senKendall, >= 10 years)
% with its p corrected for autocorrelated years (E.durAutocorr 'hamedRao': mkVarianceFactor, mkPCorrected), the
% E.durCI % interval of the slope with the same variance (senInterval), Benjamini-Hochberg over the 6 x 10 corrected p,
% and the lag-1 autocorrelation of the residuals from the Theil-Sen line (residualR1). T.p0: the classical p, T.vf: the
% variance factor (1 with 'none'). DU.timing: the same for the first and last day (DU.first, DU.last) of the three
% seasons, [10 x 3 x 2] (regions x a, c, e x onset, end), with a Benjamini-Hochberg family of its own
assert(isscalar(E.durPersist) && E.durPersist >= 1 && E.durPersist == round(E.durPersist), ...
    'E.durPersist: whole days, >= 1');
assert(isscalar(E.durMaxGap) && E.durMaxGap >= 0 && E.durMaxGap == round(E.durMaxGap), ...
    'E.durMaxGap: whole days, >= 0');
assert(numel(E.durPivot) == 2 && all(E.durPivot == round(E.durPivot)) && E.durPivot(1) >= 1 && ...
    E.durPivot(2) <= 365 && E.durPivot(1) < E.durPivot(2), 'E.durPivot: [first last] day of the 365-day calendar');
assert(isscalar(E.durAlbedo) && E.durAlbedo > 0 && E.durAlbedo < 1 && isscalar(E.durNdsi) && E.durNdsi > 0 && ...
    E.durNdsi < 1 && isscalar(E.durT2m) && isfinite(E.durT2m), ...
    'E.durAlbedo and E.durNdsi: thresholds in (0, 1); E.durT2m: a temperature in degC');
assert(isscalar(E.durPeakAlbedo) && E.durPeakAlbedo > 0 && isscalar(E.durPeakNdsi) && E.durPeakNdsi > 0 && ...
    isscalar(E.durPeakT2m) && E.durPeakT2m > 0, 'E.durPeakAlbedo, E.durPeakNdsi and E.durPeakT2m: tolerances > 0');
assert(E.durCI > 0 && E.durCI < 100, 'E.durCI: a percentage in (0, 100)');
assert(ischar(E.durAutocorr) && any(strcmp(E.durAutocorr, {'hamedRao', 'none'})), ...
    'E.durAutocorr: ''hamedRao'' or ''none''');
yb = R.years(R.ib);  nb = numel(yb);
assert(isequal(size(R.cycYear), [10 365 nb]) && isequal(X.yb(:)', yb(:)'), ...
    'stage 2g: R.cycYear and eda_extra.mat must hold the same %d baseline years', nb);
vN = find(strcmp(X.drv, 'ndsi'));  vT = find(strcmp(X.drv, 't2m'));
assert(isscalar(vN) && isscalar(vT) && size(X.drvCycYear, 3) == nb, ...
    'stage 2g: the NDSI and t2m curves of the baseline years are expected in eda_extra.mat');
assert(strcmp(X.drvUnit{vT}, 'deg C') && strcmp(X.drvUnit{vN}, '0-1'), ...
    'stage 2g: t2m in deg C and NDSI snow cover 0-1 expected');
DU.yb = yb;  DU.names = [G.regNames(G.order), {'All Alaska'}];
DU.ind = {'Albedo', 'NDSI', 't2m'};                     % the three curves
DU.metKey = {'albedo_season', 'albedo_peak', 'ndsi_season', 'ndsi_peak', 't2m_season', 't2m_peak'};
DU.metName = {sprintf('Albedo below %g', E.durAlbedo), ...
    sprintf('Albedo at its summer low (below %g and within %g of the lowest)', E.durAlbedo, E.durPeakAlbedo), ...
    sprintf('NDSI snow cover below %g', E.durNdsi), ...
    sprintf('NDSI snow cover at its summer low (below %g and within %g of the lowest)', E.durNdsi, E.durPeakNdsi), ...
    sprintf('air temperature above %g degC', E.durT2m), ...
    sprintf('air temperature at its summer high (within %g degC of the highest)', E.durPeakT2m)};
DU.metInd = [1 1 2 2 3 3];                              % the curve of each duration
below = [true true false];                              % in the state: below the threshold (Albedo, NDSI), above (t2m)
thrS = [E.durAlbedo, E.durNdsi, E.durT2m];              % the seasons' thresholds
tol = [E.durPeakAlbedo, E.durPeakNdsi, E.durPeakT2m];   % at peak: closer to the extreme than this
pw = E.durPivot(1):E.durPivot(2);
[DU.extreme, DU.pivot] = deal(nan(10, nb, 3));
[DU.thr, DU.first, DU.last, DU.val] = deal(nan(10, nb, 6));
for j = 1:10
    for i = 1:nb
        cv = {double(R.cycYear(j, :, i)), double(X.drvCycYear(j, :, i, vN)), double(X.drvCycYear(j, :, i, vT))};
        for k = 1:3
            x = cv{k};
            if mean(isfinite(x(pw))) < E.metricMinFrac, continue; end            % mid-season not observed
            xf = fillGaps(x, E.durMaxGap);
            if below(k), [xp, q] = min(xf(pw)); else, [xp, q] = max(xf(pw)); end
            p = pw(q);
            DU.extreme(j, i, k) = xp;  DU.pivot(j, i, k) = p;
            if below(k), thrP = min(thrS(k), xp + tol(k)); else, thrP = max(thrS(k), xp - tol(k)); end
            m = 2 * k - 1;                                                         % the season
            [s, e, nS] = seasonRun(xf, thrS(k), below(k), p, E.durPersist);
            DU.thr(j, i, m) = thrS(k);  DU.first(j, i, m) = s;  DU.last(j, i, m) = e;
            if isfinite(s) && isfinite(e)
                DU.val(j, i, m) = e - s + 1;
            elseif nS == 0
                DU.val(j, i, m) = 0;                                              % the pivot is out of the state
            end
            m = 2 * k;                                                             % at peak
            DU.thr(j, i, m) = thrP;
            if nS == 0
                DU.val(j, i, m) = 0;                                              % no season, so no day at its peak
            else
                [DU.val(j, i, m), DU.first(j, i, m), DU.last(j, i, m)] = stateDays(xf, thrP, below(k));
            end
        end
    end
end
% ---------- the yearly values summed up, and their trends ----------
nM = numel(DU.metKey);
[DU.med, DU.vmin, DU.vmax] = deal(nan(10, nM));          % median, lowest and highest days over the years with a value
for j = 1:10
    for m = 1:nM
        v = DU.val(j, :, m);  v = v(isfinite(v));
        if ~isempty(v), DU.med(j, m) = median(v);  DU.vmin(j, m) = min(v);  DU.vmax(j, m) = max(v); end
    end
end
T = struct('slope', nan(10, nM), 'lo', nan(10, nM), 'hi', nan(10, nM), 'tau', nan(10, nM), 'p', nan(10, nM), ...
    'q', nan(10, nM), 'n', zeros(10, nM), 'r1', nan(10, nM), 'a', nan(10, nM), 'p0', nan(10, nM), 'vf', nan(10, nM));
for j = 1:10
    for m = 1:nM
        y = DU.val(j, :, m);
        [T.slope(j, m), T.tau(j, m), T.p0(j, m), T.n(j, m), S, v] = senKendall(yb, y, 10);
        if ~isfinite(T.slope(j, m)), continue; end
        if strcmp(E.durAutocorr, 'hamedRao'), T.vf(j, m) = mkVarianceFactor(yb, y, T.slope(j, m)); else, T.vf(j, m) = 1; end
        T.p(j, m) = mkPCorrected(S, v, T.vf(j, m));
        [T.lo(j, m), T.hi(j, m)] = senInterval(yb, y, 10, E.durCI, T.vf(j, m));
        [T.r1(j, m), T.a(j, m)] = residualR1(yb, y, T.slope(j, m));
    end
end
T.q = bhFdr(T.p);                                        % Benjamini-Hochberg over the 6 durations x 10 regions
DU.trend = T;
% ---------- onset and end of the three seasons: regions x seasons (a, c, e) x [onset end] ----------
z = nan(10, 3, 2);
TM = struct('slope', z, 'lo', z, 'hi', z, 'tau', z, 'p0', z, 'vf', z, 'p', z, 'q', z, 'n', zeros(10, 3, 2), 'a', z, ...
    'med', z);
for j = 1:10
    for k = 1:3
        m = 2 * k - 1;
        for w = 1:2
            if w == 1, y = DU.first(j, :, m); else, y = DU.last(j, :, m); end
            v = y(isfinite(y));  if ~isempty(v), TM.med(j, k, w) = median(v); end
            [TM.slope(j, k, w), TM.tau(j, k, w), TM.p0(j, k, w), TM.n(j, k, w), S, v] = senKendall(yb, y, 10);
            if ~isfinite(TM.slope(j, k, w)), continue; end
            if strcmp(E.durAutocorr, 'hamedRao'), TM.vf(j, k, w) = mkVarianceFactor(yb, y, TM.slope(j, k, w));
            else, TM.vf(j, k, w) = 1; end
            TM.p(j, k, w) = mkPCorrected(S, v, TM.vf(j, k, w));
            [TM.lo(j, k, w), TM.hi(j, k, w)] = senInterval(yb, y, 10, E.durCI, TM.vf(j, k, w));
            [~, TM.a(j, k, w)] = residualR1(yb, y, TM.slope(j, k, w));
        end
    end
end
TM.q = bhFdr(TM.p);                                      % one family: onset and end of the three seasons
DU.timing = TM;
DU.settings = struct('durAlbedo', E.durAlbedo, 'durPeakAlbedo', E.durPeakAlbedo, 'durNdsi', E.durNdsi, ...
    'durPeakNdsi', E.durPeakNdsi, 'durT2m', E.durT2m, 'durPeakT2m', E.durPeakT2m, 'durPersist', E.durPersist, ...
    'durMaxGap', E.durMaxGap, 'durPivot', E.durPivot, 'durCI', E.durCI, 'durAutocorr', E.durAutocorr, 'metricMinFrac', E.metricMinFrac);
end

function [n, f1, f2] = stateDays(x, thr, below)
% the days of the curve x (1 x 365, gaps filled; NaN = no value) in the state x < thr (below = true) or x > thr (below =
% false): n their number, f1 and f2 the first and last of them. A day in the state next to a day without a value, or
% on the first or last day of the year, leaves n unknown (NaN), since the days not seen may be in the state too; such a
% first or last day leaves f1 or f2 unknown as well. No day in the state: n = 0, f1 = f2 = NaN
ok = isfinite(x);
if below, in = ok & x < thr; else, in = ok & x > thr; end
n = nnz(in);  f1 = NaN;  f2 = NaN;
if n == 0, return; end
gapB = [true, ~ok(1:end - 1)];  gapA = [~ok(2:end), true];          % no value the day before / the day after
f1 = find(in, 1, 'first');  f2 = find(in, 1, 'last');
if gapB(f1), f1 = NaN; end
if gapA(f2), f2 = NaN; end
if any(in & (gapB | gapA)), n = NaN; end
end

function y = fillGaps(x, maxGap)
% x (1 x n) with every run of <= maxGap missing values between two values filled linearly; longer runs, and runs at
% either end of x, stay NaN
y = x;  f = find(isfinite(x));
for k = find(diff(f) > 1)
    a = f(k);  b = f(k + 1);
    if b - a - 1 <= maxGap
        y(a + 1:b - 1) = x(a) + (x(b) - x(a)) * ((a + 1:b - 1) - a) / (b - a);
    end
end
end

function [s, e, nS, inS] = seasonRun(x, thr, below, p, nPersist)
% the season around day p of the curve x (1 x 365, gaps filled; NaN = no value): summer state x < thr (below = true)
% or x > thr (below = false). From the day after the last spell of nPersist consecutive winter-state days before p to
% the day before the first such spell after p (edgeSearch); shorter spells inside belong to the season. Day p not in
% the summer state: s = e = NaN and nS = 0 (no season). A boundary that cannot be found: NaN, and nS = NaN.
% nS: the summer-state days in s..e; inS: every summer-state day (logical 1 x 365)
ok = isfinite(x);
assert(ok(p), 'seasonRun: day %d has no value', p);
if below, inS = ok & x < thr; else, inS = ok & x > thr; end
st = zeros(1, numel(x));  st(inS) = 1;  st(ok & ~inS) = -1;                 % 1 summer, -1 winter, 0 no value
s = NaN;  e = NaN;  nS = 0;
if st(p) < 0, return; end
e = edgeSearch(st, p, 1, nPersist);
s = edgeSearch(st, p, -1, nPersist);
if isfinite(s) && isfinite(e), nS = nnz(inS(s:e)); else, nS = NaN; end
end

function d = edgeSearch(st, p, step, nPersist)
% from day p (summer state) in direction step (+1 later, -1 earlier): the last summer-state day before the first run of
% nPersist consecutive winter-state days (st = -1); NaN when a day without a value (st = 0) or the end of the year
% comes first
d = NaN;  last = p;  run = 0;  k = p + step;
while k >= 1 && k <= numel(st)
    if st(k) == 0, return; end
    if st(k) < 0
        run = run + 1;
        if run == nPersist, d = last; return; end
    else
        run = 0;  last = k;
    end
    k = k + step;
end
end

function [lo, hi] = senInterval(t, y, minN, level, vf)
% level % confidence interval of the Theil-Sen slope per decade (Sen 1968; Gilbert 1987, eq. 16.6, ranks rounded as in
% R's trend::sens.slope): the N pairwise slopes sorted, C = z((1 + level / 100) / 2) x sqrt(Var S) with the Mann-Kendall
% variance of senKendall (ties in y corrected) times vf (optional, default 1: the Hamed-Rao factor of
% mkVarianceFactor); the slopes of rank round((N - C) / 2) and round((N + C) / 2) + 1. NaN with < minN
% values or when a rank falls outside 1..N
if nargin < 5 || ~isfinite(vf), vf = 1; end
ok = isfinite(y(:));  t = t(:);  t = t(ok);  y = y(:);  y = y(ok);  n = numel(y);
lo = NaN;  hi = NaN;
if n < minN, return; end
[i, j] = find(triu(true(n), 1));
d = sort(10 * (y(j) - y(i)) ./ (t(j) - t(i)));  N = numel(d);
[~, ~, g] = unique(y);  tp = accumarray(g, 1);
v = (n * (n - 1) * (2 * n + 5) - sum(tp .* (tp - 1) .* (2 * tp + 5))) / 18;
C = sqrt(2) * erfinv(level / 100) * sqrt(v * vf);
r1 = round((N - C) / 2);  r2 = round((N + C) / 2) + 1;
if r1 >= 1 && r2 <= N, lo = d(r1);  hi = d(r2); end
end

function [r, a] = residualR1(t, y, slopeDec)
% lag-1 autocorrelation of the residuals from the Theil-Sen line (slope slopeDec per decade; intercept a = median of
% y - slope x t over the years with a value, Conover) over the pairs of consecutive years that both have a value
% (residuals centred on their mean); NaN without a slope or with < 3 such pairs
ok = isfinite(y(:));  t = t(:);  y = y(:);  r = NaN;  a = NaN;
if ~isfinite(slopeDec), return; end
b = slopeDec / 10;  a = median(y(ok) - b * t(ok));
d = y - a - b * t;  d(ok) = d(ok) - mean(d(ok));
pr = ok(1:end - 1) & ok(2:end) & diff(t) == 1;
if nnz(pr) < 3, return; end
r = sum(d([pr; false]) .* d([false; pr])) / sum(d(ok) .^ 2);
end

function f = mkVarianceFactor(t, y, slopeDec)
% Hamed and Rao (1998) factor n/n* on the Mann-Kendall variance for autocorrelated years: the residuals
% from the Theil-Sen line (slope slopeDec per decade) over the n years with a value, ranked (ties: mean ranks) and
% centred; rho_k = sum over the pairs of years k apart that both have a value / the sum of squares, k = 1 .. span;
% f = 1 + 2 / (n (n - 1) (n - 2)) x sum of (n - k)(n - k - 1)(n - k - 2) rho_k over the k with |rho_k| >
% 1.96 / sqrt(n) (with no missing year this equals pymannkendall's hamed_rao_modification_test). Kept >= 1: with
% ~23 years the rank autocorrelations are noisy, and a factor < 1 would make the test more liberal than the classical
% one. NaN without a slope or < 3 values
ok = isfinite(y(:));  t = t(:);  t = t(ok);  y = y(:);  y = y(ok);  n = numel(y);  f = NaN;
if ~isfinite(slopeDec) || n < 3, return; end
d = y - slopeDec / 10 * (t - mean(t));
d = round(d * 1e6) / 1e6;                         % equal residuals tie exactly (the median slope joins two years)
[ds, o] = sort(d);  rk = zeros(n, 1);  i = 1;
while i <= n                                       % mean ranks of ties
    k = i;  while k < n && ds(k + 1) == ds(i), k = k + 1; end
    rk(o(i:k)) = (i + k) / 2;  i = k + 1;
end
rk = rk - mean(rk);  den = sum(rk .^ 2);  f = 1;
if den == 0, return; end
zc = sqrt(2) * erfinv(0.95);  s = 0;               % 1.96
for k = 1:n - 1
    [hit, w] = ismember(t + k, t);
    if ~any(hit), continue; end
    rho = sum(rk(hit) .* rk(w(hit))) / den;
    if abs(rho) > zc / sqrt(n), s = s + (n - k) * (n - k - 1) * (n - k - 2) * rho; end
end
f = max(1, 1 + 2 / (n * (n - 1) * (n - 2)) * s);
end

function p = mkPCorrected(S, v, f)
% two-sided Mann-Kendall p with the variance v of the score S multiplied by f (continuity corrected as senKendall)
p = NaN;
if ~isfinite(S) || ~isfinite(v) || ~isfinite(f) || v <= 0, return; end
z = (S - sign(S)) / sqrt(v * f);  if S == 0, z = 0; end
p = erfc(abs(z) / sqrt(2));
end

function printDurations(DU, E)
nb = numel(DU.yb);  T = DU.trend;  nM = numel(DU.metKey);
fprintf(['Durations %d-%d (stage 2g; spells >= %d days end a season, gaps <= %d days filled): median days over ' ...
    'the years with a value (years), Theil-Sen days per decade (* Benjamini-Hochberg q < 0.05, %d tests; ' ...
    'Mann-Kendall p: %s)\n'], DU.yb(1), DU.yb(end), E.durPersist, E.durMaxGap, nnz(isfinite(T.p)), ...
    strrep(strrep(DU.settings.durAutocorr, 'hamedRao', 'Hamed-Rao corrected'), 'none', 'classical'));
lab = {sprintf('Albedo < %g', E.durAlbedo), 'Albedo at low', sprintf('NDSI < %g', E.durNdsi), 'NDSI at low', ...
    sprintf('t2m > %g', E.durT2m), 't2m at high'};
fprintf('  %-30s', 'region');  fprintf('%-18s', lab{:});  fprintf('\n');
for j = [10 1:9]
    fprintf('  %-30s', DU.names{j});
    for m = 1:nM
        if T.n(j, m) == 0, md = '-'; else, md = sprintf('%g', DU.med(j, m)); end
        if isfinite(T.slope(j, m)), sl = sprintf('%+.1f', T.slope(j, m)); else, sl = '-'; end
        if isfinite(T.q(j, m)) && T.q(j, m) < 0.05, sl = [sl '*']; end
        fprintf('%-18s', sprintf('%s (%d) %s', md, T.n(j, m), sl));
    end
    fprintf('\n');
end
TM = DU.timing;  sn = {'Albedo season', 'NDSI snow-free season', 'thaw season'};
fprintf(['  Onset and end (days per decade, + later; Benjamini-Hochberg over their %d tests, %d with q < 0.05; ' ...
    'smallest q %.2f)\n'], nnz(isfinite(TM.p)), nnz(TM.q < 0.05), min(TM.q(:)));
for k = 1:3
    fprintf('    %-22s All Alaska onset %+5.1f, end %+5.1f', sn{k}, TM.slope(10, k, 1), TM.slope(10, k, 2));
    qk = TM.q(:, k, :);
    if any(isfinite(qk(:)))
        [qm, i] = min(qk(:));  [jj, ~, ww] = ind2sub([10 1 2], i);  wn = {'onset', 'end'};
        fprintf('; smallest q %.2f: %s %s %+.1f', qm, DU.names{jj}, wn{ww}, TM.slope(jj, k, ww));
    end
    fprintf('\n');
end
for m = 1:nM
    k = DU.metInd(m);  hasCurve = isfinite(DU.pivot(:, :, k));
    unk = hasCurve & isnan(DU.val(:, :, m));
    fprintf(['  %s: days unknown in %d of %d region-years (first day unknown in %d, last day in %d), no summer ' ...
        'curve in %d\n'], DU.metName{m}, nnz(unk), 10 * nb, nnz(unk & isnan(DU.first(:, :, m))), ...
        nnz(unk & isnan(DU.last(:, :, m))), nnz(~hasCurve));
end
end

function drawDurationFigures(DU, G, E)
% fig18a-fig18f (E.figures 18; one figure per duration), All Alaska across the top, then the nine ecoregions; the panels
% of a figure share their axis limits, so the regions compare directly. Days at peak (b, d, f): the days in every year
% with the Theil-Sen line (solid: Benjamini-Hochberg q < 0.05 over the tests of the six durations; dashed: not
% significant). Seasons (a, c, e): each year a bar from the onset (circle) to the end (square) on a day-of-year axis
% (its length is the duration), the Theil-Sen lines of the onset and the end (solid: q < 0.05 in DU.timing's family),
% the three trends in the subtitles
if isempty(DU) || ~any(E.figures == 18), return; end
C = palette();  names = [G.short(G.order), {'All Alaska'}];  P = [10 1:9];  yb = DU.yb;  T = DU.trend;
S = DU.settings;
cInd = [C.ink; 0.106 0.604 0.667; 0.851 0.345 0.169];         % Albedo, NDSI, t2m (the colours of fig17)
xt = yb(mod(yb, 5) == 0);  xl = [yb(1) - 0.5, yb(end) + 0.5];
dg = '^{\circ}C';                                             % layoutTitle: TeX
head = {sprintf('Albedo below %g: days from its spring to its autumn crossing', S.durAlbedo), ...
    sprintf('Albedo at its summer low: days within %g of the year''s lowest Albedo (and below %g)', ...
        S.durPeakAlbedo, S.durAlbedo), ...
    sprintf('NDSI snow cover below %g (snow-free season): days from its spring to its autumn crossing', S.durNdsi), ...
    sprintf('NDSI snow cover at its summer low: days within %g of the year''s lowest NDSI (and below %g)', ...
        S.durPeakNdsi, S.durNdsi), ...
    sprintf('Air temperature above %g%s (thaw season): days from its spring to its autumn crossing', S.durT2m, dg), ...
    sprintf('Air temperature at its summer high: days within %g%s of the year''s highest (the warmest weeks)', ...
        S.durPeakT2m, dg)};
rule = {sprintf(['%d-%d, 7-day smoothed regional curves; the season ends when the curve stays across the ' ...
    'threshold for %d or more days.'], yb(1), yb(end), S.durPersist), ...
    sprintf('%d-%d, 7-day smoothed regional curves; every day within the tolerance counts.', yb(1), yb(end))};
if strcmp(S.durAutocorr, 'hamedRao'), mkNote = ', Mann-Kendall p and CI corrected for autocorrelation (Hamed-Rao)';
else, mkNote = ''; end
trendKey = sprintf(['Theil-Sen trends (days per decade): solid if Benjamini-Hochberg q < 0.05 across %d tests%s; ' ...
    'dashed otherwise.'], nnz(isfinite(T.p)), mkNote);
TM = DU.timing;
timeKey = sprintf(['Bars: each year from its onset (circle) to its end (square), so the bar length is the days. ' ...
    'Lines: Theil-Sen trends of the onset and the end, solid if Benjamini-Hochberg q < 0.05 across their %d tests, ' ...
    'dashed otherwise.\nSubtitles: the trend of the days (q across the %d duration tests), then of the onset and the ' ...
    'end, in days per decade (+ = later)%s.'], nnz(isfinite(TM.p)), nnz(isfinite(T.p)), ...
    strrep(mkNote, ', Mann-Kendall', '; Mann-Kendall'));
cOn = [0.106 0.604 0.667];  cEnd = [0.851 0.345 0.169];  cBar = [0.78 0.78 0.78];   % onset, end, the bars
ms = [1 32 60 91 121 152 182 213 244 274 305 335 366];                             % the first of each month (365 days)
mn = {'1 Jan', '1 Feb', '1 Mar', '1 Apr', '1 May', '1 Jun', '1 Jul', '1 Aug', '1 Sep', '1 Oct', '1 Nov', '1 Dec', ...
    '31 Dec'};
ndsiNote = {sprintf(['No MODIS NDSI from early October: where NDSI is still below %g then, the end is unknown ' ...
    '(subtitles: years with a value).'], S.durNdsi), ['No MODIS NDSI from early October: where NDSI is still at ' ...
    'its low then, the days are unknown (subtitles: years with a value).']};
for m = 1:numel(DU.metKey)
    col = cInd(DU.metInd(m), :);
    v = DU.val(:, :, m);  v = v(isfinite(v));
    if isempty(v) || max(v) == 0, yl = [0 10]; else, yl = [0, 10 * ceil(1.02 * max(v) / 10)]; end
    kind = 2 - mod(m, 2);                                             % 1 a season, 2 at peak
    if kind == 1                                                      % a season: the day-of-year axis
        k = (m + 1) / 2;
        d = [DU.first(:, :, m), DU.last(:, :, m)];  d = d(isfinite(d));
        if isempty(d), yl = [121 305]; else, yl = [ms(find(ms <= min(d), 1, 'last')), ms(find(ms >= max(d), 1))]; end
        yt = ms(ms >= yl(1) & ms <= yl(2));  ytl = mn(ms >= yl(1) & ms <= yl(2));
        str = sprintf('%s\n%s\n%s', head{m}, rule{kind}, timeKey);
    else
        yt = 0:50:200;  ytl = arrayfun(@num2str, yt, 'UniformOutput', false);
        str = sprintf('%s\n%s\n%s', head{m}, rule{kind}, trendKey);
    end
    if DU.metInd(m) == 2, str = sprintf('%s\n%s', str, ndsiNote{kind}); end
    f = newFig(13, 10.5, E);
    t = tiles(f, 4, 3);
    layoutTitle(t, str, C);
    for p = 1:10
        j = P(p);
        if p == 1, a = tileAt(t, 1, [1 3], C); else, a = tileAt(t, p + 2, [1 1], C); end
        if kind == 1
            f1 = DU.first(j, :, m);  f2 = DU.last(j, :, m);  ok = isfinite(f1) & isfinite(f2);
            if any(ok)
                xb = [yb(ok); yb(ok); nan(1, nnz(ok))];  yv = [f1(ok); f2(ok); nan(1, nnz(ok))];
                plot(a, xb(:), yv(:), '-', 'Color', cBar, 'LineWidth', 4, 'Marker', 'none');
            end
            plot(a, yb, f1, 'o', 'Color', cOn, 'LineStyle', 'none', 'MarkerSize', 4, 'MarkerFaceColor', cOn);
            plot(a, yb, f2, 's', 'Color', cEnd, 'LineStyle', 'none', 'MarkerSize', 4, 'MarkerFaceColor', cEnd);
            cw = [cOn; cEnd];
            for w = 1:2
                if ~isfinite(TM.slope(j, k, w)), continue; end
                if TM.q(j, k, w) < 0.05, st = '-'; else, st = '--'; end
                x = [yb(1) yb(end)];
                plot(a, x, TM.a(j, k, w) + TM.slope(j, k, w) / 10 * x, 'Color', cw(w, :), 'LineStyle', st, ...
                    'LineWidth', 1.2, 'Marker', 'none');
            end
        else
            plot(a, yb, DU.val(j, :, m), '-o', 'Color', col, 'LineWidth', 0.8, 'MarkerSize', 4, 'MarkerFaceColor', col);
        end
        if kind == 2 && isfinite(T.slope(j, m))
            if T.q(j, m) < 0.05, st = '-'; else, st = '--'; end
            x = [yb(1) yb(end)];
            plot(a, x, T.a(j, m) + T.slope(j, m) / 10 * x, 'Color', col, 'LineStyle', st, 'LineWidth', 1.2, ...
                'Marker', 'none');
        end
        set(a, 'XLim', xl, 'YLim', yl, 'XGrid', 'on', 'YGrid', 'on');

        if p == 1 || p > 7
            set(a, 'XTick', xt)
        else
            set(a, 'XTick', xt, 'XTickLabel', []);
        end
        if p == 1 || p == 2 || p == 5 || p == 8
            set(a, 'YTick', yt, 'YTickLabel', ytl)
        else
            set(a, 'YTick', yt, 'YTickLabel', []);
        end
        if p==1
            set(a, 'XTick', 2003:2025);
        end
        panelTitle(a, names{j}, C, durPanelText(DU, j, m));
    end
    if kind == 1
        ylabel(t, 'Day of year', 'FontSize',13,'Interpreter','latex');
    else
        ylabel(t, 'Duration (days)', 'FontSize',13,'Interpreter','latex');
    end
    saveFig(f, fullfile(E.figDir, sprintf('fig18%s_%s.pdf', char('a' + m - 1), DU.metKey{m})));
end
end

function s = durPanelText(DU, j, m)
% the subtitle lines of region j in the figure of duration m: the Theil-Sen slope with its E.durCI % interval,
% 'Sen +2.1 d/dec [95% CI +0.8, +3.4]' ('Duration ...' for a season), and the Benjamini-Hochberg q with the years that
% have a value, 'q = 0.79; n = 23/23 years'; or why there is no trend. A season adds the trends of its onset and end
% (timingText). The median is in tab18b
T = DU.trend;  nb = numel(DU.yb);  n = T.n(j, m);
nStr = sprintf('n = %d/%d years', n, nb);
if ~isfinite(T.slope(j, m))
    s = {'no trend (fewer than 10 years)', nStr};
    if mod(m, 2) == 1, s{3} = timingText(DU, j, m); end
    return
end
if mod(m, 2) == 1, lab = 'Duration'; else, lab = 'Sen'; end           % a season: three trends
if isfinite(T.lo(j, m)) && isfinite(T.hi(j, m))
    s1 = sprintf('%s %+.1f d/dec [%g%% CI %+.1f, %+.1f]', lab, T.slope(j, m), DU.settings.durCI, T.lo(j, m), ...
        T.hi(j, m));
else
    s1 = sprintf('%s %+.1f d/dec [CI not available]', lab, T.slope(j, m));
end
s = {s1, sprintf('%s; %s', qText(T.q(j, m)), nStr)};
if mod(m, 2) == 1, s{3} = timingText(DU, j, m); end
end

function s = qText(q)
% 'q = 0.79', 'q < 0.001' or 'q not available' (subtitles of fig18a-f)
if ~isfinite(q), s = 'q not available'; elseif q < 0.001, s = 'q < 0.001'; else, s = sprintf('q = %.2g', q); end
end

function s = timingText(DU, j, m)
% third subtitle line of a season (m = 1, 3, 5) of region j: 'Onset -7.9 d/dec (q = 0.67), end +14.6 d/dec
% (q = 0.18)'; a part without a trend (fewer than 10 years) says so: 'Onset +0.0 d/dec (q = 1), end: < 10 years'
TM = DU.timing;  k = (m + 1) / 2;  wn = {'Onset', 'end'};  p = cell(1, 2);  ok = isfinite(squeeze(TM.slope(j, k, :)));
for w = 1:2
    if ok(w)
        p{w} = sprintf('%s %+.1f d/dec (%s)', wn{w}, TM.slope(j, k, w), qText(TM.q(j, k, w)));
    else
        p{w} = sprintf('%s: < 10 years', wn{w});
    end
end
s = sprintf('%s, %s', p{1}, p{2});
end

function writeDurationTables(DU, G, E)
% tab18: every region and year, each curve's summer extreme and its day, and each duration's threshold, first and last
% day and days; tab18b: per region and duration the years with a value, the median, lowest and highest days and the
% trend. Every run with E.durations
if isempty(DU), return; end
full = [G.regNames(G.order), {'All Alaska'}];  yb = DU.yb;  nb = numel(yb);  nM = numel(DU.metKey);
extName = {'albedo_summer_low', 'ndsi_summer_low', 't2m_summer_high'};  un = {'', '', '_degC'};
hdr = {'region', 'year'};
for k = 1:3
    hdr = [hdr, {[extName{k} un{k}], [extName{k} '_doy']}];                                           %#ok<AGROW>
    for m = find(DU.metInd == k)
        hdr = [hdr, {[DU.metKey{m} '_threshold' un{k}], [DU.metKey{m} '_first_doy'], [DU.metKey{m} '_last_doy'], ...
            [DU.metKey{m} '_days']}];                                                                    %#ok<AGROW>
    end
end
rows = cell(10 * nb, numel(hdr));  q = 0;
for j = 1:10
    for i = 1:nb
        r = {full{j}, yb(i)};
        for k = 1:3
            r = [r, {DU.extreme(j, i, k), DU.pivot(j, i, k)}];                                          %#ok<AGROW>
            for m = find(DU.metInd == k)
                r = [r, {DU.thr(j, i, m), DU.first(j, i, m), DU.last(j, i, m), DU.val(j, i, m)}];       %#ok<AGROW>
            end
        end
        q = q + 1;  rows(q, :) = r;
    end
end
writeCsv(fullfile(E.tabDir, 'tab18_season_lengths.csv'), hdr, rows);
T = DU.trend;  rows = cell(10 * nM, 16);  q = 0;
for j = 1:10
    for m = 1:nM
        q = q + 1;
        rows(q, :) = {full{j}, DU.metKey{m}, DU.metName{m}, T.n(j, m), DU.med(j, m), DU.vmin(j, m), DU.vmax(j, m), ...
            T.slope(j, m), T.lo(j, m), T.hi(j, m), T.tau(j, m), T.p0(j, m), T.vf(j, m), T.p(j, m), T.q(j, m), ...
            T.r1(j, m)};
    end
end
writeCsv(fullfile(E.tabDir, 'tab18b_season_trends.csv'), {'region', 'duration', 'definition', 'years', ...
    'median_days', 'min_days', 'max_days', 'theil_sen_days_per_decade', sprintf('ci%g_low', E.durCI), ...
    sprintf('ci%g_high', E.durCI), 'kendall_tau', 'mann_kendall_p_classical', 'mk_variance_factor', ...
    'mann_kendall_p', 'bh_fdr_q', 'lag1_autocorr_residuals'}, rows);
TM = DU.timing;  rows = cell(10 * 3 * 2, 14);  q = 0;  wn = {'onset', 'end'};
for j = 1:10
    for k = 1:3
        for w = 1:2
            q = q + 1;
            rows(q, :) = {full{j}, DU.metKey{2 * k - 1}, wn{w}, TM.n(j, k, w), TM.med(j, k, w), TM.slope(j, k, w), ...
                TM.lo(j, k, w), TM.hi(j, k, w), TM.tau(j, k, w), TM.p0(j, k, w), TM.vf(j, k, w), TM.p(j, k, w), ...
                TM.q(j, k, w), TM.a(j, k, w)};
        end
    end
end
writeCsv(fullfile(E.tabDir, 'tab18c_season_onset_end_trends.csv'), {'region', 'season', 'crossing', 'years', ...
    'median_doy', 'theil_sen_days_per_decade', sprintf('ci%g_low', E.durCI), sprintf('ci%g_high', E.durCI), ...
    'kendall_tau', 'mann_kendall_p_classical', 'mk_variance_factor', 'mann_kendall_p', 'bh_fdr_q_timing', ...
    'theil_sen_intercept_doy'}, rows);
fprintf('  wrote tab18, tab18b and tab18c to %s\n', E.tabDir);
end

% ------------------- RQ3: temperature sensitivity and the temperature-scaled estimates of E.whatIfYears (stage 2h) -------------------

function T = readWhatIfCsv(fp, cols, isText)
% CMIP6 CSV table (cmip6_assemble.py) -> struct with a field per column (text: cell; numbers: double, empty -> NaN).
% The header must hold exactly the columns cols, in this order; every row must parse
assert(isfile(fp), 'CMIP6 table %s not found (written by cmip6_assemble.py; README.md in E.cmip6Dir)', fp);
txt = strrep(fileread(fp), sprintf('\r'), '');
if ~isempty(txt) && txt(end) ~= sprintf('\n'), txt = [txt sprintf('\n')]; end
nl = find(txt == sprintf('\n'));
assert(numel(nl) >= 2, 'CMIP6 table %s holds no data rows', fp);
hdr = regexprep(strsplit(txt(1:nl(1) - 1), ','), '^"(.*)"$', '$1');   % names with or without quotes
assert(isequal(hdr, cols), 'CMIP6 table %s has the columns %s; expected %s', fp, strjoin(hdr, ','), strjoin(cols, ','));
fmt = repmat({'%f'}, 1, numel(cols));  fmt(isText) = {'%s'};
body = txt(nl(1) + 1:end);
nRows = numel(nl) - 1;
Cc = textscan(body, strjoin(fmt, ''), 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
for c = 1:numel(cols)
    assert(numel(Cc{c}) == nRows, 'CMIP6 table %s: column %s has %d values for %d rows', fp, cols{c}, ...
        numel(Cc{c}), nRows);
    v = Cc{c};
    if isText(c), v = regexprep(v, '^"(.*)"$', '$1'); end
    T.(cols{c}) = v;
end
end

function s = scenLabel(sc)
% 'ssp245' -> 'SSP2-4.5'
s = sprintf('SSP%s-%s.%s', sc(4), sc(5), sc(6));
end

function fn = cmip6Table(E, s, kind)
% file name of a CMIP6 table of stage 2h for scenario E.whatIfScenarios{s} (README.md of E.cmip6Dir):
% kind 'period' the block means of the baseline and of E.whatIfYears, alb_cmip6_tas_<scenario>_<first>-<last>.csv;
% 'years' the block means of every year, alb_cmip6_tas_years_<scenario>.csv
switch kind
    case 'period'
        fn = sprintf('alb_cmip6_tas_%s_%d-%d.csv', E.whatIfScenarios{s}, E.whatIfYears(1), E.whatIfYears(2));
    case 'years'
        fn = sprintf('alb_cmip6_tas_years_%s.csv', E.whatIfScenarios{s});
    otherwise
        error('cmip6Table: kind ''period'' or ''years'', not ''%s''', kind);
end
end

function AT = edaAlbTemp(R, G, E, X)
% Stage 2h (RQ3): how strongly the Albedo of the sunlit season responds to April-May air temperature from year to year
% in every region (rows as R.cycMean: the 9 ecoregions in G.order, then All Alaska), and a what-if for E.whatIfYears
% that applies the April-May warming projected by NEX-GDDP-CMIP6 to that response (a temperature-scaled estimate for
% every year of E.whatIfYears and for their mean). From each baseline year's 7-day smoothed regional curves (those of
% fig01: Albedo R.cycYear, 2-m air temperature X.drvCycYear; gaps of <= E.albTempMaxGap days filled linearly, fillGaps)
% and the daily climatology of the incoming shortwave (R.swIn.mean, ERA5-Land SSRD, MJ m-2 d-1):
%   season   the days of the 8-day blocks E.albTempBlocks (365-day calendar)
%   A        the season's sunlight-weighted Albedo: the sum over its days of w(d) x Albedo(d), w(d) = the SSRD
%            climatology of day d / its sum over the season (the same weights in every year), so 1 - A = the share of
%            today's mean incoming shortwave of the season that the surface absorbs; NaN unless every day has a value
%   T        the mean 2-m air temperature over the days of the blocks E.albTempTBlocks; NaN unless every day has a value
%   fit      A = a + b T over the years with both (>= E.albTempMinYears): b the sensitivity (Albedo per K), its OLS
%            standard error, r, the two-sided p of b (t with n - 2 degrees of freedom), b with the year as a second
%            covariate (linear trends removed), the lag-1 autocorrelation of the residuals of consecutive years, and the
%            E.albTempCI % interval of b: percentiles over E.albTempNBoot resamples of the years (rng E.albTempSeed,
%            twister; the same resampled years for every region; a resample refits the line on its years, a repeated
%            year counted again). Also: the E.albTempCI % t interval of b (n - 2 degrees of freedom); the range of b
%            with each year left out; the year's coefficient with T held fixed (per decade, its OLS standard error and
%            p, n - 3 degrees of freedom); b with the season's backup-inversion share or valid share (albTempQA) as a
%            second covariate (with its OLS standard error); and the parts of b from the blocks before, in and after the
%            T window: each block's sunlight-weighted Albedo regressed on T, times the block's weight w_k (they add up
%            to b, as A = the sum of w_k x the blocks' Albedo)
%   blocks   check of the season sum: per 8-day block k of the season, the slope beta_k of the block's sunlight-weighted
%            Albedo on the mean temperature of the E.albTempWindow blocks ending with k; U = the sum of w_k beta_k (w_k
%            the block's share of the season's sunlight, so A = the sum of w_k x the block's Albedo): the season's
%            response to 1 K of warming in every block; in the what-if, the sum of w_k beta_k dT_k with each block's own
%            projected warming, which lets every block (also those of autumn) respond to its own temperature
%   what-if  per scenario of E.whatIfScenarios, model and region: dT = the warming of the blocks E.albTempTBlocks (the
%            mean over E.whatIfYears minus that over the baseline: albTempCmip6), dA = b dT, and (delta method) T + dT
%            and A + dA from the means of the fitted years; the extra absorbed shortwave dQ = -dA x the season's sum of
%            the SSRD climatology (MJ m-2 per year), also as a mean flux over the whole year (W m-2). Central value: the
%            median of the models with b; ranges: E.albTempCI % of every model with every resample (pooled), of the
%            models alone (with b) and of the resamples alone (with the models' mean warming). Also the mean of the
%            models: dTMean, dAMean = b x their mean warming, aFutMean, dQMean and dQWMean (the values of fig19c and
%            tab19c; the medians stay in tab19)
%   yearly   per scenario, model, region and year y of E.whatIfYears: dT(y) = the mean over the blocks
%            E.albTempTBlocks of the block's mean in year y minus its baseline mean (albTempCmip6), and A(y) = the
%            mean A of the fitted years + b dT(y), the temperature-scaled Albedo of year y (the same b every year).
%            Per year over the models: the mean, median, quartiles and E.albTempCI % range of A(y), the mean, median
%            and range of dT(y), and how many models are warmer than the warmest fitted year (fig19, tab19d; every
%            model in tab19e). The mean of the yearly means is the period's aFutMean, as each block of T holds the
%            same days every year
% Fields of AT: yb, names, blocks, tBlocks, window, days and tDays ([first last] day of the season and of T), scen,
% years, yf (the years of E.whatIfYears, first to last), pixels [9 x 2] (codes 1-9: pixels of the grid, pixels with a
% model value), draws [E.albTempNBoot x years] (uint8: the resampled years), settings, qaShare and qaValid [10 x years]
% (the season's backup-inversion share and valid share), and
%   Q, qYear, share [10 x 1]  the SSRD climatology summed over the season and over the year (MJ m-2), their ratio
%   alb, t [10 x years]  A and T;  albBlk, tWin [10 x blocks x years]  each block's Albedo and window temperature;
%       wBlk [10 x blocks]  the blocks' shares of the season's sunlight
%   fit.<f> [10 x 1]  n, slope, icp (a), se, r, p, lo, hi, slopeDet, r1, meanA, sdA, meanT, sdT, minT, maxT, minA,
%       maxA, tLo, tHi, looMin, looMax, yearTrend, yearSE (per decade), yearP, slopeShare, slopeShareSE, slopeValid,
%       slopeValidSE; fit.part [10 x 3] (b from the blocks before, in, after the T window); fit.used [10 x years]
%   boot.slope, boot.icp [10 x resamples]  b and a of every resample
%   blk  n [10 x 1], beta [10 x blocks], U, ULo, UHi [10 x 1], betaBoot [10 x blocks x resamples]
%   S(s)  models {1 x M}, excluded, dTb [M x 10 x 46] (block warming, K; NaN outside the exported blocks), dT [M x 10],
%       dTw [M x 10 x blocks] (the warming of each block's window), dA and dAblk [M x 10] (b dT; sum of w_k beta_k
%       dT_k), sum.<f> [10 x 1] (tab19, writeAlbTempTables); dTy and aY [M x 10 x years] (the warming and the Albedo
%       A(y) of every model and year of E.whatIfYears; aY NaN without a fit) and yr.<f> [10 x years] (dTMean, dTMed,
%       dTLo, dTHi, aMean, aMed, aQ1, aQ3, aLo, aHi, dAMean, nWarmer; tab19d)
Kb = E.albTempBlocks(:)';  Kt = E.albTempTBlocks(:)';  win = E.albTempWindow;  nB = E.albTempNBoot;
assert(isscalar(win) && win >= 1 && win == round(win), 'E.albTempWindow: whole 8-day blocks, >= 1');
assert(~isempty(Kb) && isequal(Kb, Kb(1):Kb(end)) && Kb(1) >= win && Kb(end) <= 45, ['E.albTempBlocks: consecutive ' ...
    '8-day blocks from %d (the temperature window of the first block reaches %d blocks back) to 45'], win, win - 1);
assert(~isempty(Kt) && isequal(Kt, Kt(1):Kt(end)) && Kt(1) >= 1 && Kt(end) <= 45, ...
    'E.albTempTBlocks: consecutive 8-day blocks within 1-45');
assert(isscalar(E.albTempMaxGap) && E.albTempMaxGap >= 0 && E.albTempMaxGap == round(E.albTempMaxGap), ...
    'E.albTempMaxGap: whole days, >= 0');
assert(E.albTempCI > 0 && E.albTempCI < 100 && nB >= 100 && nB <= 1e5 && nB == round(nB), ...
    'E.albTempCI: a percentage in (0, 100); E.albTempNBoot: a whole number of resamples, 100 to 100,000');
yb = R.years(R.ib);  yb = yb(:)';  nb = numel(yb);
assert(E.albTempMinYears >= 3 && E.albTempMinYears <= nb, 'E.albTempMinYears: 3 to the %d baseline years', nb);
assert(isequal(size(R.cycYear), [10 365 nb]) && isequal(X.yb(:)', yb), ...
    'stage 2h: R.cycYear and eda_extra.mat must hold the same %d baseline years', nb);
vT = find(strcmp(X.drv, 't2m'));
assert(isscalar(vT) && size(X.drvCycYear, 3) == nb && strcmp(X.drvUnit{vT}, 'deg C'), ...
    'stage 2h: the 2-m air temperature curves (deg C) of the baseline years are expected in eda_extra.mat');
assert(isequal(size(R.swIn.mean), [10 365]), 'stage 2h: R.swIn.mean [10 x 365] expected in eda_results.mat');
nK = numel(Kb);  pLo = (100 - E.albTempCI) / 2;  pHi = 100 - pLo;
dS = 8 * Kb(1) - 7 : 8 * Kb(end);              % days of the season (365-day calendar)
dT = 8 * Kt(1) - 7 : 8 * Kt(end);              % days of T
kw = Kb(1) - win + 1 : Kb(end);                % blocks of the temperature windows of the season's blocks
dW = 8 * kw(1) - 7 : 8 * kw(end);
AT.yb = yb;  AT.names = [G.regNames(G.order), {'All Alaska'}];  AT.blocks = Kb;  AT.tBlocks = Kt;  AT.window = win;
AT.days = dS([1 end]);  AT.tDays = dT([1 end]);
% ---------- the season's sunlight-weighted Albedo and the temperature T of every region and year ----------
[AT.Q, AT.qYear, AT.share] = deal(nan(10, 1));
[AT.alb, AT.t] = deal(nan(10, nb));
[AT.albBlk, AT.tWin] = deal(nan(10, nK, nb));
AT.wBlk = nan(10, nK);
for j = 1:10
    sw = R.swIn.mean(j, dS);
    assert(all(isfinite(sw)) && all(sw > 0), 'stage 2h: %s: no incoming shortwave climatology on a day of the season', ...
        AT.names{j});
    w = sw / sum(sw);  wk = sum(reshape(w, 8, nK), 1);
    AT.wBlk(j, :) = wk;  AT.Q(j) = sum(sw);
    if all(isfinite(R.swIn.mean(j, :))), AT.qYear(j) = sum(R.swIn.mean(j, :));  AT.share(j) = AT.Q(j) / AT.qYear(j); end
    for i = 1:nb
        a = fillGaps(double(R.cycYear(j, :, i)), E.albTempMaxGap);
        if all(isfinite(a(dS)))
            AT.alb(j, i) = sum(w .* a(dS));
            ab = sum(reshape(w .* a(dS), 8, nK), 1) ./ wk;          % each block's sunlight-weighted Albedo
            assert(abs(sum(wk .* ab) - AT.alb(j, i)) <= 1e-12, 'stage 2h: the blocks'' Albedo does not add up to A');
            AT.albBlk(j, :, i) = ab;
        end
        x = fillGaps(double(X.drvCycYear(j, :, i, vT)), E.albTempMaxGap);
        if all(isfinite(x(dT))), AT.t(j, i) = mean(x(dT)); end
        tb = mean(reshape(x(dW), 8, numel(kw)), 1);                    % block means (NaN unless every day)
        AT.tWin(j, :, i) = arrayfun(@(q) mean(tb(q:q + win - 1)), 1:nK);   % the window ending with block Kb(q)
    end
end
[AT.qaShare, AT.qaValid] = albTempQA(G, E, yb, dS);   % the season's retrieval mix
% ---------- the fit, its resamples and the block-by-block check ----------
rng(E.albTempSeed, 'twister');
DR = randi(nb, nB, nb);                          % row d: the years of resample d (indices into yb)
AT.draws = uint8(DR);
Cn = zeros(nB, nb);                              % how often each year is drawn in each resample
for i = 1:nb, Cn(:, i) = sum(DR == i, 2); end
F = struct();
for f = {'slope', 'icp', 'se', 'r', 'p', 'lo', 'hi', 'slopeDet', 'r1', 'meanA', 'sdA', 'meanT', 'sdT', 'minT', 'maxT', ...
        'minA', 'maxA', 'tLo', 'tHi', 'looMin', 'looMax', 'yearTrend', 'yearSE', 'yearP', 'slopeShare', ...
        'slopeShareSE', 'slopeValid', 'slopeValidSE'}
    F.(f{1}) = nan(10, 1);
end
F.n = zeros(10, 1);  F.used = false(10, nb);  F.part = nan(10, 3);
AT.boot = struct('slope', nan(10, nB), 'icp', nan(10, nB));
B = struct('n', zeros(10, 1), 'beta', nan(10, nK), 'U', nan(10, 1), 'ULo', nan(10, 1), 'UHi', nan(10, 1), ...
    'betaBoot', nan(10, nK, nB));
for j = 1:10
    m = isfinite(AT.alb(j, :)) & isfinite(AT.t(j, :));  n = nnz(m);  F.n(j) = n;
    if n < E.albTempMinYears, continue; end
    ctx = sprintf('stage 2h, %s', AT.names{j});
    A = AT.alb(j, m)';  T = AT.t(j, m)';  yy = yb(m)';
    [b, a0, e] = slopeFit(T, A, ctx);
    df = n - 2;  se = sqrt(sum(e .^ 2) / df / sum((T - mean(T)) .^ 2));
    F.used(j, :) = m;  F.slope(j) = b;  F.icp(j) = a0;  F.se(j) = se;  F.r(j) = pearson(T, A);
    F.p(j) = tPvalue(b / se, df);
    q = tQuantile(E.albTempCI, df);  F.tLo(j) = b - q * se;  F.tHi(j) = b + q * se;
    bl = arrayfun(@(i) slopeFit(T([1:i - 1, i + 1:n]), A([1:i - 1, i + 1:n]), ctx), 1:n);   % each year left out
    F.looMin(j) = min(bl);  F.looMax(j) = max(bl);
    [bd, sd] = olsFit([ones(n, 1), T - mean(T), yy - mean(yy)], A - mean(A), [ctx ', with the year']);
    F.slopeDet(j) = bd(2);  F.yearTrend(j) = 10 * bd(3);  F.yearSE(j) = 10 * sd(3);
    F.yearP(j) = tPvalue(bd(3) / sd(3), n - 3);
    sh = AT.qaShare(j, m)';  va = AT.qaValid(j, m)';               % the retrieval mix as a second covariate
    [bq, sq] = olsFit([ones(n, 1), T - mean(T), sh - mean(sh)], A - mean(A), [ctx ', with the backup-inversion share']);
    F.slopeShare(j) = bq(2);  F.slopeShareSE(j) = sq(2);
    [bq, sq] = olsFit([ones(n, 1), T - mean(T), va - mean(va)], A - mean(A), [ctx ', with the valid share']);
    F.slopeValid(j) = bq(2);  F.slopeValidSE(j) = sq(2);
    pc = AT.wBlk(j, :) .* arrayfun(@(k) slopeFit(T, reshape(AT.albBlk(j, k, m), [], 1), ctx), 1:nK);   % w_k b_k
    assert(abs(sum(pc) - b) <= 1e-9 * max(abs(b), 1e-3), '%s: the blocks'' parts do not add up to the slope', ctx);
    F.part(j, :) = [sum(pc(Kb < Kt(1))), sum(pc(Kb >= Kt(1) & Kb <= Kt(end))), sum(pc(Kb > Kt(end)))];
    cons = find(diff(yy) == 1);                                      % residuals of consecutive years
    if numel(cons) >= 3, F.r1(j) = pearson(e(cons), e(cons + 1)); end
    F.meanA(j) = mean(A);  F.sdA(j) = std(A);  F.meanT(j) = mean(T);  F.sdT(j) = std(T);
    F.minT(j) = min(T);  F.maxT(j) = max(T);  F.minA(j) = min(A);  F.maxA(j) = max(A);
    c = Cn(:, m);
    assert(all(sum(c > 0, 2) >= 3), '%s: a resample holds fewer than 3 of the fitted years', ctx);
    [AT.boot.slope(j, :), AT.boot.icp(j, :)] = slopeBoot(c, T, A);
    F.lo(j) = pctl(AT.boot.slope(j, :), pLo);  F.hi(j) = pctl(AT.boot.slope(j, :), pHi);
    % block by block, over the fitted years whose window temperatures all have a value
    m2 = m & all(isfinite(reshape(AT.tWin(j, :, :), nK, nb)), 1);  B.n(j) = nnz(m2);
    if B.n(j) < E.albTempMinYears, continue; end
    c2 = Cn(:, m2);
    assert(all(sum(c2 > 0, 2) >= 3), '%s: a resample holds fewer than 3 of the years of the block check', ctx);
    for q = 1:nK
        Aq = reshape(AT.albBlk(j, q, m2), [], 1);  Tq = reshape(AT.tWin(j, q, m2), [], 1);
        B.beta(j, q) = slopeFit(Tq, Aq, sprintf('%s, block %d', ctx, Kb(q)));
        B.betaBoot(j, q, :) = reshape(slopeBoot(c2, Tq, Aq), 1, 1, nB);
    end
    B.U(j) = AT.wBlk(j, :) * B.beta(j, :)';
    ub = AT.wBlk(j, :) * reshape(B.betaBoot(j, :, :), nK, nB);
    B.ULo(j) = pctl(ub, pLo);  B.UHi(j) = pctl(ub, pHi);
end
AT.fit = F;  AT.blk = B;
% ---------- the what-if: the warming of NEX-GDDP-CMIP6 x the sensitivity ----------
[Sc, AT.pixels] = albTempCmip6(G, E, yb, unique([Kt, kw]), Kt);
AT.scen = E.whatIfScenarios;  AT.years = E.whatIfYears;  AT.yf = E.whatIfYears(1):E.whatIfYears(2);  nYf = numel(AT.yf);
kq = 1e6 / (365 * 86400);                        % MJ m-2 per year -> mean W m-2 over the year
for s = 1:numel(Sc)
    S = Sc(s);  M = numel(S.models);  ctx = sprintf('stage 2h, %s', scenLabel(AT.scen{s}));
    S.dT = mean(S.dTb(:, :, Kt), 3);                                  % [M x 10]
    assert(size(S.dTy, 1) == M && size(S.dTy, 2) == 10 && size(S.dTy, 3) == nYf && ...
        max(abs(reshape(mean(S.dTy, 3) - S.dT, 1, []))) <= 1e-4, ['%s: the mean of the yearly warming differs from ' ...
        'the warming of E.whatIfYears'], ctx);                       % the two tables agree
    S.dTw = nan(M, 10, nK);
    for q = 1:nK, S.dTw(:, :, q) = mean(S.dTb(:, :, Kb(q) - win + 1 : Kb(q)), 3); end
    assert(all(isfinite(S.dT(:))) && all(isfinite(S.dTw(:))), '%s: a model without the warming of a block', ctx);
    [S.dA, S.dAblk] = deal(nan(M, 10));
    sm = struct();
    for f = {'nModels', 'dT', 'dTLo', 'dTHi', 'dTsd', 'asWarm', 'tFut', 'tFutLo', 'tFutHi', 'dA', 'dALo', 'dAHi', ...
            'dAmLo', 'dAmHi', 'dAbLo', 'dAbHi', 'aFut', 'aFutLo', 'aFutHi', 'dQ', 'dQLo', 'dQHi', 'dQW', 'dQWLo', ...
            'dQWHi', 'dAblk', 'dAblkLo', 'dAblkHi', 'dTMean', 'dAMean', 'aFutMean', 'dQMean', 'dQWMean'}
        sm.(f{1}) = nan(10, 1);
    end
    for j = 1:10
        d = S.dT(:, j);                                              % [M x 1] the warming of every model (K)
        sm.nModels(j) = M;  sm.dT(j) = pctl(d, 50);  sm.dTLo(j) = pctl(d, pLo);  sm.dTHi(j) = pctl(d, pHi);
        sm.dTMean(j) = mean(d);
        if ~isfinite(F.slope(j)), continue; end
        sm.dTsd(j) = sm.dT(j) / F.sdT(j);
        sm.asWarm(j) = nnz(AT.t(j, F.used(j, :)) >= F.meanT(j) + sm.dT(j));
        sm.tFut(j) = F.meanT(j) + sm.dT(j);  sm.tFutLo(j) = F.meanT(j) + sm.dTLo(j);  sm.tFutHi(j) = F.meanT(j) + sm.dTHi(j);
        S.dA(:, j) = F.slope(j) * d;
        X1 = d * AT.boot.slope(j, :);                                % [M x resamples]: every model with every resample
        sm.dA(j) = pctl(S.dA(:, j), 50);  sm.dALo(j) = pctl(X1(:), pLo);  sm.dAHi(j) = pctl(X1(:), pHi);
        sm.dAmLo(j) = pctl(S.dA(:, j), pLo);  sm.dAmHi(j) = pctl(S.dA(:, j), pHi);
        xb = mean(d) * AT.boot.slope(j, :);                          % the resamples with the models' mean warming
        sm.dAbLo(j) = pctl(xb, pLo);  sm.dAbHi(j) = pctl(xb, pHi);
        sm.aFut(j) = F.meanA(j) + sm.dA(j);  sm.aFutLo(j) = F.meanA(j) + sm.dALo(j);  sm.aFutHi(j) = F.meanA(j) + sm.dAHi(j);
        sm.dAMean(j) = F.slope(j) * mean(d);  sm.aFutMean(j) = F.meanA(j) + sm.dAMean(j);   % the mean of the models
        sm.dQMean(j) = -sm.dAMean(j) * AT.Q(j);  sm.dQWMean(j) = sm.dQMean(j) * kq;
        q1 = -S.dA(:, j) * AT.Q(j);  qp = -X1(:) * AT.Q(j);           % extra absorbed shortwave (MJ m-2 per year)
        sm.dQ(j) = pctl(q1, 50);  sm.dQLo(j) = pctl(qp, pLo);  sm.dQHi(j) = pctl(qp, pHi);
        sm.dQW(j) = sm.dQ(j) * kq;  sm.dQWLo(j) = sm.dQLo(j) * kq;  sm.dQWHi(j) = sm.dQHi(j) * kq;
        if isfinite(B.U(j))
            dw = reshape(S.dTw(:, j, :), M, nK);                     % [M x blocks]
            S.dAblk(:, j) = dw * (AT.wBlk(j, :) .* B.beta(j, :))';
            Xb = (reshape(B.betaBoot(j, :, :), nK, nB) .* AT.wBlk(j, :)')' * dw';   % [resamples x M]
            sm.dAblk(j) = pctl(S.dAblk(:, j), 50);  sm.dAblkLo(j) = pctl(Xb(:), pLo);  sm.dAblkHi(j) = pctl(Xb(:), pHi);
        end
    end
    S.sum = sm;
    % ---------- every year of E.whatIfYears: A(y) = the mean A + b dT(y), every model ----------
    S.aY = nan(M, 10, nYf);  yr = struct();
    for f = {'dTMean', 'dTMed', 'dTLo', 'dTHi', 'aMean', 'aMed', 'aQ1', 'aQ3', 'aLo', 'aHi', 'dAMean', 'nWarmer'}
        yr.(f{1}) = nan(10, nYf);
    end
    for j = 1:10
        for i = 1:nYf
            d = S.dTy(:, j, i);                                      % [M x 1] the warming of year i (K)
            yr.dTMean(j, i) = mean(d);  yr.dTMed(j, i) = pctl(d, 50);  yr.dTLo(j, i) = pctl(d, pLo);
            yr.dTHi(j, i) = pctl(d, pHi);
            if ~isfinite(F.slope(j)), continue; end
            a = F.meanA(j) + F.slope(j) * d;  S.aY(:, j, i) = a;
            yr.aMean(j, i) = mean(a);  yr.aMed(j, i) = pctl(a, 50);  yr.aQ1(j, i) = pctl(a, 25);
            yr.aQ3(j, i) = pctl(a, 75);  yr.aLo(j, i) = pctl(a, pLo);  yr.aHi(j, i) = pctl(a, pHi);
            yr.dAMean(j, i) = F.slope(j) * mean(d);
            yr.nWarmer(j, i) = nnz(F.meanT(j) + d > F.maxT(j));    % warmer than the warmest fitted year
        end
    end
    S.yr = yr;
    AT.S(s) = S;
end
AT.settings = struct('albTempBlocks', Kb, 'albTempTBlocks', Kt, 'albTempWindow', win, 'albTempMaxGap', ...
    E.albTempMaxGap, 'albTempMinYears', E.albTempMinYears, 'albTempNBoot', nB, 'albTempSeed', E.albTempSeed, ...
    'albTempCI', E.albTempCI, 'whatIfScenarios', {E.whatIfScenarios}, 'whatIfYears', E.whatIfYears, ...
    'whatIfMinDays', E.whatIfMinDays);
end

function [b, a, e] = slopeFit(x, y, what)
% least-squares line y = a + b x of the column vectors x and y (what: the fit, named in the error message); e: the
% residuals
xm = mean(x);  ym = mean(y);  dx = x - xm;  sxx = sum(dx .^ 2);
assert(sxx > 0, '%s: the temperature is the same in every year', what);
b = sum(dx .* (y - ym)) / sxx;  a = ym - b * xm;  e = y - a - b * x;
end

function [b, a] = slopeBoot(c, x, y)
% least-squares line of every resample of the years: c [resamples x years] holds how often each year (x, y: column
% vectors) is drawn, so a resample's line is the fit weighted by these counts (a repeated year counted again);
% b, a: [1 x resamples]. Centred first, so the sums do not cancel
x0 = mean(x);  y0 = mean(y);  dx = x - x0;  dy = y - y0;
s = sum(c, 2);  mx = (c * dx) ./ s;  my = (c * dy) ./ s;
sxx = c * (dx .^ 2) - s .* mx .^ 2;  sxy = c * (dx .* dy) - s .* mx .* my;
assert(all(sxx > 0), 'slopeBoot: the years of a resample all have the same temperature');
b = (sxy ./ sxx)';  a = (y0 + my - (sxy ./ sxx) .* (x0 + mx))';
end

function p = tPvalue(t, df)
% two-sided p of Student's t with df degrees of freedom: I(df / (df + t^2); df/2, 1/2)
p = betainc(df / (df + t ^ 2), df / 2, 0.5);
end

function [beta, se] = olsFit(Xd, y, what)
% least squares (lsq) with the classical standard errors sqrt(diag(s2 (X'X)^-1)), s2 = the residual sum of squares /
% (n - columns); what: the fit, named in error messages
beta = lsq(Xd, y, what);
[n, k] = size(Xd);
assert(n > k, '%s: no degrees of freedom left', what);
e = y - Xd * beta;  s2 = sum(e .^ 2) / (n - k);
[~, Rq] = qr(Xd, 0);  Ri = Rq \ eye(k);                  % (X'X)^-1 = R^-1 R^-T
se = sqrt(s2 * sum(Ri .^ 2, 2));
end

function [sh, va] = albTempQA(G, E, yb, dS)
% stage 2h: the retrieval mix of the season (days dS of the 365-day calendar) in every region (G.order, then All
% Alaska = codes 1-10, the rows of R.cycYear) and baseline year yb, from the region x day sums of stage 1 (alb_YYYY.mat
% in E.cacheDir; 29 Feb dropped as in to365): sh = the season's backup-inversion share of the valid values (the sum of
% rmag, each value's share of its area from MCD43 magnitude inversions 1 - QA/100, over the sum of rcnt), va = the
% season's valid share of the pixel-days (the sum of rcnt over pixels x days) [10 x years]
nb = numel(yb);  [sh, va] = deal(nan(10, nb));
npx = [reshape(G.npix(G.order + 1), 1, []), sum(G.npix(2:11))];
for i = 1:nb
    fp = fullfile(E.cacheDir, sprintf('alb_%d.mat', yb(i)));
    L = load(fp, 'rcnt', 'rmag');
    assert(isfield(L, 'rcnt') && isfield(L, 'rmag') && isequal(size(L.rcnt), [11 366]) && ...
        isequal(size(L.rmag), [11 366]), 'stage 2h: %s lacks rcnt and rmag [11 x 366]', fp);
    assert(all(L.rmag(:) >= 0 & L.rmag(:) <= L.rcnt(:) * (1 + 1e-12)), ['stage 2h: %s: more backup inversions than ' ...
        'values on a day (rmag > rcnt) or a negative rmag'], fp);
    c = to365([L.rcnt(G.order + 1, :); sum(L.rcnt(2:11, :), 1)], yb(i));
    g = to365([L.rmag(G.order + 1, :); sum(L.rmag(2:11, :), 1)], yb(i));
    c = sum(c(:, dS), 2);  g = sum(g(:, dS), 2);
    assert(all(c > 0), 'stage 2h: %s: a region without a valid value in the season', fp);
    sh(:, i) = g ./ c;  va(:, i) = c ./ (npx(:) * numel(dS));
end
end

function [S, pix] = albTempCmip6(G, E, yb, need, kt)
% stage 2h: the CMIP6 warming per scenario (E.whatIfScenarios), model, region (G.order, then All Alaska) and 8-day
% block, read from the tables in E.cmip6Dir (README.md there): every model, ecoregion code and exported block once; a
% model is used if every exported block keeps >= E.whatIfMinDays of its days in both periods; All Alaska = the mean of
% the 9 ecoregions weighted by their pixels with a model value. need: the blocks stage 2h uses (they must be in the
% tables); kt: the blocks of T (E.albTempTBlocks). S(s): models, excluded, dTb [M x 10 x 46] (mean over E.whatIfYears
% minus mean over the baseline years yb, K; NaN outside the exported blocks), dTy [M x 10 x years] (the warming of every
% year of E.whatIfYears, the mean over the blocks kt of the block's mean in that year minus its baseline mean, K); pix
% [9 x 2]: pixels of codes 1-9 in the grid and with a model value. Two tables per scenario (cmip6Table): the period
% table (block means of the baseline and of E.whatIfYears) and the yearly table (the block means of every year). For the
% models used, the yearly table must hold every ecoregion code, block of kt and year of the baseline and of
% E.whatIfYears once, each block with the same days in every year of E.whatIfYears, and its day-weighted means over the
% baseline and over E.whatIfYears must equal the period table (1e-4 K, the same days): the two tables must come from the
% same runs
scen = E.whatIfScenarios;
assert(iscell(scen) && ~isempty(scen) && all(cellfun(@(c) ~isempty(regexp(c, '^ssp\d\d\d$', 'once')), scen)), ...
    'E.whatIfScenarios: a cell of scenario names such as ''ssp245''');
assert(numel(E.whatIfYears) == 2 && E.whatIfYears(2) >= E.whatIfYears(1) && E.whatIfYears(1) > yb(end), ...
    'E.whatIfYears: [first last] year, after the baseline');
assert(E.whatIfMinDays > 0 && E.whatIfMinDays <= 1, 'E.whatIfMinDays in (0, 1]');
kt = kt(:)';
assert(~isempty(kt) && all(ismember(kt, need)), 'stage 2h: the blocks of T must be among the blocks it reads');
nYb = numel(yb);  nYf = E.whatIfYears(2) - E.whatIfYears(1) + 1;  nk = numel(kt);
yAll = [yb(:)', E.whatIfYears(1):E.whatIfYears(2)];  iB = 1:nYb;  iF = nYb + 1 : nYb + nYf;   % years of the yearly table
fp = fullfile(E.cmip6Dir, 'alb_cmip6_pixels.csv');
P = readWhatIfCsv(fp, {'code', 'pixels', 'pixels_with_tas'}, false(1, 3));
assert(isequal(sort(P.code(:))', 1:9), 'stage 2h: %s must hold the ecoregion codes 1-9, once each', fp);
[~, o] = sort(P.code);  pix = [P.pixels(o), P.pixels_with_tas(o)];
npix = reshape(G.npix(2:10), 1, []);                        % G.npix(c + 1): pixels of code c
assert(isequal(pix(:, 1)', npix), 'stage 2h: the pixels per ecoregion in %s (%s) differ from the analysis grid (%s)', ...
    fp, num2str(pix(:, 1)'), num2str(npix));
assert(all(pix(:, 2) <= pix(:, 1)) && all(pix(:, 2) > 0), ['stage 2h: %s: every ecoregion needs pixels with a model ' ...
    'value, no more than its pixels'], fp);
wc = reshape(pix(:, 2), 1, 9);                              % All Alaska: weights of codes 1-9
cols = {'model', 'scenario', 'code', 'block', 'tas_base_K', 'tas_future_K', 'n_base_days', 'n_future_days'};
colsY = {'model', 'scenario', 'code', 'block', 'year', 'tas_K', 'n_days'};
S = struct('models', {}, 'excluded', {}, 'dTb', {}, 'dTy', {});
for s = 1:numel(scen)
    fp = fullfile(E.cmip6Dir, cmip6Table(E, s, 'period'));
    C = readWhatIfCsv(fp, cols, [true true false(1, 6)]);
    assert(all(strcmp(C.scenario, scen{s})), 'stage 2h: %s holds rows of another scenario', fp);
    assert(all(ismember(C.code, 1:9)) && all(ismember(C.block, 1:46)), 'stage 2h: %s: codes 1-9 and blocks 1-46', fp);
    assert(isequal(unique(C.code(:))', 1:9), 'stage 2h: %s lacks ecoregion code(s) %s', fp, ...
        num2str(setdiff(1:9, C.code(:)')));
    kb = unique(C.block(:))';                               % the exported blocks
    assert(all(ismember(need, kb)), ['stage 2h: %s lacks block(s) %s; stage 2h needs blocks %d-%d (the CMIP6 ' ...
        'tables hold blocks 5-37: README.md in E.cmip6Dir)'], fp, num2str(setdiff(need, kb)), min(need), max(need));
    models = unique(C.model);  models = models(:)';  M = numel(models);
    [~, im] = ismember(C.model, models);
    sub = [im(:), C.code(:), C.block(:)];
    cnt = accumarray(sub, 1, [M 9 46]);
    assert(all(reshape(cnt(:, :, kb), 1, []) == 1) && ~any(reshape(cnt(:, :, setdiff(1:46, kb)), 1, [])), ...
        'stage 2h: %s must hold every model, ecoregion code and block exactly once', fp);
    put = @(v) accumarray(sub, v(:), [M 9 46], [], NaN);
    Tb = put(C.tas_base_K);  Tf = put(C.tas_future_K);  Nb = put(C.n_base_days);  Nf = put(C.n_future_days);
    assert(max(reshape(Nb(:, :, kb), 1, [])) == 8 * nYb && max(reshape(Nf(:, :, kb), 1, [])) == 8 * nYf, ['stage 2h: %s: ' ...
        'no block holds %d baseline and %d future days (8 days x %d and %d years): the table was made for other ' ...
        'periods than the baseline (%d-%d) and E.whatIfYears'], fp, 8 * nYb, 8 * nYf, nYb, nYf, yb(1), yb(end));
    use = true(1, M);
    for i = 1:M
        nbd = reshape(Nb(i, :, kb), 1, []);  nfd = reshape(Nf(i, :, kb), 1, []);
        assert(all(nbd <= 8 * nYb) && all(nfd <= 8 * nYf), ['stage 2h: %s, %s: more days in a block than %d and %d ' ...
            'years hold'], fp, models{i}, nYb, nYf);
        use(i) = all(nbd >= E.whatIfMinDays * 8 * nYb) && all(nfd >= E.whatIfMinDays * 8 * nYf);
    end
    assert(nnz(use) >= 3, 'stage 2h: %s: fewer than 3 models with >= %g %% of the days of every block', fp, ...
        100 * E.whatIfMinDays);
    Tb = Tb(use, :, :);  Tf = Tf(use, :, :);  Nb = Nb(use, :, :);  Nf = Nf(use, :, :);
    v = [reshape(Tb(:, :, kb), 1, []), reshape(Tf(:, :, kb), 1, [])];
    assert(all(isfinite(v) & v > 190 & v < 320), 'stage 2h: %s: block mean temperatures outside 190-320 K', fp);
    d = Tf - Tb;                                            % [M x 9 codes x 46]
    dTb = nan(nnz(use), 10, 46);
    dTb(:, 1:9, :) = d(:, G.order, :);
    dTb(:, 10, :) = sum(d .* wc, 2) / sum(wc);
    % ---------- every year of E.whatIfYears, from the yearly table ----------
    fy = fullfile(E.cmip6Dir, cmip6Table(E, s, 'years'));
    Y = readWhatIfCsv(fy, colsY, [true true false(1, 5)]);
    assert(all(strcmp(Y.scenario, scen{s})), 'stage 2h: %s holds rows of another scenario', fy);
    mu = models(use);  Mu = numel(mu);
    [inM, im] = ismember(Y.model, mu);  [inK, ik] = ismember(Y.block, kt);  [inY, iy] = ismember(Y.year, yAll);
    r = inM(:) & inK(:) & inY(:);                           % rows of the models used, blocks kt and years needed
    assert(all(ismember(Y.code(r), 1:9)), 'stage 2h: %s: ecoregion codes 1-9 expected', fy);
    sub = [im(r), Y.code(r), ik(r), iy(r)];
    cnt = accumarray(sub, 1, [Mu 9 nk numel(yAll)]);
    assert(all(cnt(:) == 1), ['stage 2h: %s must hold every model used, ecoregion code, block %d-%d and year (%d-%d ' ...
        'and %d-%d) exactly once (README.md of E.cmip6Dir)'], fy, kt(1), kt(end), yb(1), yb(end), E.whatIfYears(1), ...
        E.whatIfYears(2));
    Ty = accumarray(sub, Y.tas_K(r), size(cnt));  Ny = accumarray(sub, Y.n_days(r), size(cnt));
    assert(all(isfinite(Ty(:)) & Ty(:) > 190 & Ty(:) < 320), 'stage 2h: %s: yearly block means outside 190-320 K', fy);
    assert(all(Ny(:) >= 1 & Ny(:) <= 8 & Ny(:) == round(Ny(:))), 'stage 2h: %s: 1-8 days per block and year expected', fy);
    assert(all(reshape(Ny(:, :, :, iF) == Ny(:, :, :, iF(1)), 1, [])), ['stage 2h: %s: a block holds different ' ...
        'numbers of days in different years of E.whatIfYears'], fy);
    PT = {Tb(:, :, kt), Nb(:, :, kt), iB, 'the baseline years';  Tf(:, :, kt), Nf(:, :, kt), iF, 'E.whatIfYears'};
    for q = 1:2                                             % the yearly table against the period table
        n = sum(Ny(:, :, :, PT{q, 3}), 4);  m = sum(Ty(:, :, :, PT{q, 3}) .* Ny(:, :, :, PT{q, 3}), 4) ./ n;
        assert(isequal(n, PT{q, 2}) && max(abs(m(:) - reshape(PT{q, 1}, [], 1))) <= 1e-4, ['stage 2h: the day-' ...
            'weighted mean of %s over %s differs from %s by > 1e-4 K, or holds other days: the two tables do not ' ...
            'come from the same runs (README.md of E.cmip6Dir)'], fy, PT{q, 4}, fp);
    end
    dy = reshape(mean(Ty(:, :, :, iF) - Tb(:, :, kt), 3), Mu, 9, nYf);   % [M x 9 codes x years] (K)
    dTy = nan(Mu, 10, nYf);
    dTy(:, 1:9, :) = dy(:, G.order, :);
    dTy(:, 10, :) = sum(dy .* wc, 2) / sum(wc);
    S(s) = struct('models', {mu}, 'excluded', {models(~use)}, 'dTb', dTb, 'dTy', dTy);
end
end

function printAlbTemp(AT, E)
F = AT.fit;  B = AT.blk;
fprintf(['RQ3 (stage 2h): sunlight-weighted Albedo of %s-%s against the mean air temperature of %s-%s, %d-%d (%g %% ' ...
    'intervals from %d resamples of the years; t: Student''s t); temperature-scaled estimate for %d-%d (NEX-GDDP-CMIP6 ' ...
    'warming x slope; medians and means of the models, ranges of models x resamples):\n'], doyLabel(AT.days(1)), ...
    doyLabel(AT.days(2)), doyLabel(AT.tDays(1)), doyLabel(AT.tDays(2)), AT.yb(1), AT.yb(end), E.albTempCI, ...
    E.albTempNBoot, AT.years(1), AT.years(2));
for j = [10 1:9]
    if ~isfinite(F.slope(j))
        fprintf('  %s: %d years with both, fewer than %d: no fit\n', AT.names{j}, F.n(j), E.albTempMinYears);  continue
    end
    fprintf(['  %s (%d years; %.0f %% of the year''s sunlight): Albedo %.3f (SD %.3f), T %+.1f degC (SD %.1f K); slope ' ...
        '%+.4f per K (%+.4f to %+.4f; t %+.4f to %+.4f; each year left out %+.4f to %+.4f), r %+.2f, p %.1g; ' ...
        'residual lag-1 r %+.2f\n'], AT.names{j}, F.n(j), 100 * AT.share(j), F.meanA(j), F.sdA(j), F.meanT(j), ...
        F.sdT(j), F.slope(j), F.lo(j), F.hi(j), F.tLo(j), F.tHi(j), F.looMin(j), F.looMax(j), F.r(j), F.p(j), F.r1(j));
    fprintf(['    checks: with the year %+.4f per K, the year %+.4f per decade (SE %.4f, p %.2g); with the backup-' ...
        'inversion share %+.4f (SE %.4f), with the valid share %+.4f (SE %.4f); parts before / in / after the T ' ...
        'window %+.4f / %+.4f / %+.4f; block by block %+.4f per K (%+.4f to %+.4f)\n'], F.slopeDet(j), ...
        F.yearTrend(j), F.yearSE(j), F.yearP(j), F.slopeShare(j), F.slopeShareSE(j), F.slopeValid(j), ...
        F.slopeValidSE(j), F.part(j, 1), F.part(j, 2), F.part(j, 3), B.U(j), B.ULo(j), B.UHi(j));
    for s = 1:numel(AT.S)
        m = AT.S(s).sum;
        fprintf(['    %s (%d models): warming %+.1f K (%+.1f to %+.1f; %.1f SD; %d of %d years as warm); Albedo %+.4f ' ...
            '(%+.4f to %+.4f), %.3f in %d-%d; mean of the models %+.1f K, %+.4f, %.3f; extra absorbed %+.0f MJ m-2 a ' ...
            'year = %+.2f W m-2 (%+.2f to %+.2f); block by block %+.4f (%+.4f to %+.4f)\n'], scenLabel(AT.scen{s}), ...
            m.nModels(j), m.dT(j), m.dTLo(j), m.dTHi(j), m.dTsd(j), m.asWarm(j), F.n(j), m.dA(j), m.dALo(j), ...
            m.dAHi(j), m.aFut(j), AT.years(1), AT.years(2), m.dTMean(j), m.dAMean(j), m.aFutMean(j), m.dQ(j), ...
            m.dQW(j), m.dQWLo(j), m.dQWHi(j), m.dAblk(j), m.dAblkLo(j), m.dAblkHi(j));
    end
end
if isfinite(F.slope(10))                                    % every year: All Alaska
    for s = 1:numel(AT.S)
        y = AT.S(s).yr;  nm = numel(AT.S(s).models);
        fprintf(['  %s, All Alaska, every year: Albedo, mean of the models (model-years warmer than the warmest fitted ' ...
            'year: %d of %d)\n    %s\n'], scenLabel(AT.scen{s}), sum(y.nWarmer(10, :)), nm * numel(AT.yf), ...
            strjoin(arrayfun(@(i) sprintf('%d %.3f', AT.yf(i), y.aMean(10, i)), 1:numel(AT.yf), ...
            'UniformOutput', false), ', '));
    end
end
for s = 1:numel(AT.S)
    if ~isempty(AT.S(s).excluded)
        fprintf('  %s: models left out (a block with < %g %% of its days): %s\n', scenLabel(AT.scen{s}), ...
            100 * E.whatIfMinDays, strjoin(AT.S(s).excluded, ', '));
    end
end
end

function drawAlbTempFigure(AT, G, E)
% fig19 (E.figures 19; RQ3): All Alaska, one panel per scenario (E.whatIfScenarios, top to bottom) on one time axis,
% from the first baseline year - 1 to the last year of E.whatIfYears + 1, and one Albedo axis:
%   observed   the season's sunlight-weighted Albedo of every fitted baseline year (dots joined by a line) and its mean
%              (dashed, across the whole panel; its value in the key)
%   estimate   the years of E.whatIfYears on a shaded background: for every year the temperature-scaled Albedo of
%              every model (A(y) = the mean + b dT(y); tab19d, tab19e): a box over the middle half of the models
%              (25-75 %), whiskers over the E.albTempCI % range, a tick at their median, and a line with dots through
%              the yearly means of the models, all in the scenario's colour (IPCC AR6 WGI colours of the SSPs:
%              scenColour; the box filled with a light tint of it, the line darker)
%   subtitle   the number of models and the E.whatIfYears mean of the models with its change from the baseline mean
%              (aFutMean and dAMean of tab19: the mean of the yearly means)
% One key at the bottom, in neutral grey (the colour names the scenario, the panel title too). The grid is drawn under
% the data ('Layer' 'bottom': on top, MATLAB draws grid lines across the boxes). The boxes show how the models'
% April-May warming differs from model to model and from year to year (each model year is one realisation of the
% weather); the uncertainty of b is not in them (its E.albTempCI % range is in tab19; the period ranges of models x
% resamples in tab19 and tab19c)
if isempty(AT) || ~any(E.figures == 19), return; end
C = palette();  nS = numel(AT.scen);  F = AT.fit;  j = 10;  yb = AT.yb;  yf = AT.yf;  nYf = numel(yf);
assert(isfinite(F.slope(j)), 'fig19: no fit for All Alaska');
pLo = (100 - E.albTempCI) / 2;  pHi = 100 - pLo;
sl = cellfun(@scenLabel, AT.scen, 'UniformOutput', false);  col = scenColour(AT.scen);
by = sprintf('%d-%d', yb(1), yb(end));  fy = sprintf('%d-%d', yf(1), yf(end));
ss = sprintf('%s - %s', doyLabel(AT.days(1)), doyLabel(AT.days(2)));
mth = @(fmt, x) ['$', sprintf(fmt, x), '$'];                       % numbers in math mode: true minus signs
used = F.used(j, :);  y = yb(used);  A = AT.alb(j, used);  a0 = F.meanA(j);   % the fitted years
Y = arrayfun(@(s) AT.S(s).yr, 1:nS);
for s = 1:nS
    v = [Y(s).aMean(j, :), Y(s).aMed(j, :), Y(s).aQ1(j, :), Y(s).aQ3(j, :), Y(s).aLo(j, :), Y(s).aHi(j, :)];
    assert(numel(v) == 6 * nYf && all(isfinite(v)), 'fig19: %s has no estimate for a year of All Alaska', sl{s});
end
yl = niceLimits([A, a0, arrayfun(@(s) min(Y(s).aLo(j, :)), 1:nS), arrayfun(@(s) max(Y(s).aHi(j, :)), 1:nS)]);
xl = [yb(1) - 1, yf(end) + 1];  wb = 0.6;                          % years; box width (years)
f = newFig(11, 1.2 + 2.3 * nS, E);
t = tiles(f, nS, 1);
for s = 1:nS
    a = tileAt(t, s, [1 1], C);
    q = Y(s);  c = col(s, :);  cf = 0.25 * c + 0.75;  cd = 0.75 * c;   % the scenario's colour, its tint and dark shade
    patch(a, [yf(1) - 0.5, yf(end) + 0.5, yf(end) + 0.5, yf(1) - 0.5], yl([1 1 2 2]), [0.5961    0.7843    0.9373], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.2, 'HandleVisibility','off');                                         % the years of the estimate, shaded
    ho = plot(a, y, A, '-o', 'Color', C.ink, 'LineWidth', 1, 'MarkerSize', 3.5, 'MarkerFaceColor', 'w', 'DisplayName',sprintf('Observed, %s', by));
    hm = plot(a, xl, a0 * [1 1], '--', 'Color', C.ink2, 'LineWidth', 0.8, 'DisplayName',sprintf('Mean observed (%s)', mth('%.3f', a0)));
    for i = 1:nYf
        x = yf(i);
        plot(a, [x x], [q.aLo(j, i), q.aQ1(j, i)], '-', 'Color', C.ink, 'LineWidth', 0.8, 'HandleVisibility','off');     % whiskers
        plot(a, [x x], [q.aQ3(j, i), q.aHi(j, i)], '-', 'Color', C.ink, 'LineWidth', 0.8, 'HandleVisibility','off');
        patch(a, x + wb / 2 * [-1 1 1 -1], [q.aQ1(j, i), q.aQ1(j, i), q.aQ3(j, i), q.aQ3(j, i)], 'w', 'EdgeColor', C.ink, ...
            'LineWidth', 0.7, 'HandleVisibility','off');                                     % the middle half of the models
        plot(a, x + wb / 2 * [-1 1], q.aMed(j, i) * [1 1], '-', 'Color', C.ink, 'LineWidth', 1.4, 'HandleVisibility','off');   % their median
    end
    plot(a, yf, q.aMean(j, :), '-o', 'Color', [0.0588    0.4510    0.9608], 'LineWidth', 1, 'MarkerSize', 3.5, 'MarkerFaceColor', [0.0588    0.4510    0.9608], 'MarkerEdgeColor', 'none' ,'DisplayName', 'Mean of the models');
    m = AT.S(s).sum;
    panelTitle(a, sl{s}, C, sprintf('%d models; %s mean of the models: %s (%s from the %s mean)', ...
        numel(AT.S(s).models), fy, mth('%.3f', m.aFutMean(j)), mth('%+.3f', m.dAMean(j)), by));
    set(a, 'XLim', xl, 'YLim', yl, 'XTick', 5 * ceil(xl(1) / 5) : 5 : xl(2), 'XGrid', 'on', 'YGrid', 'on', ...
        'Layer', 'bottom', 'TickLabelInterpreter','latex');
    if s < nS, set(a, 'XTickLabel', {}); end %else, xlabel(a, '', 'Interpreter', 'latex'); 
end
ylabel(t, 'Blue Sky Albedo', 'Interpreter', 'latex');
g = C.ink2;  gf = 0.25 * g + 0.75;                                 % the key: marks in neutral grey
% hk = [ho, hm, plot(a, NaN, NaN, '-o', 'Color', g, 'LineWidth', 1.6, 'MarkerSize', 3.5, 'MarkerFaceColor', g), ...
%       patch(a, NaN(1, 4), NaN(1, 4), gf, 'EdgeColor', g, 'LineWidth', 0.7), ...
%       plot(a, NaN, NaN, '-', 'Color', g, 'LineWidth', 0.8), plot(a, NaN, NaN, '-', 'Color', g, 'LineWidth', 1.4)];
% lab = {sprintf('observed, %s', by), sprintf('mean %s (%s)', by, mth('%.3f', a0)), ...
%     'mean of the models (colour: the scenario)', 'box: middle half of the models (25-75 %)', ...
%     sprintf('whiskers: %g %% of the models (%g-%g %%)', E.albTempCI, pLo, pHi), 'tick: median of the models'};
lg = legend(a, 'show');
set(lg, 'FontSize', 8.5, 'Interpreter', 'latex', 'Orientation', 'horizontal', ...
    'NumColumns', 3, 'Box', 'on');
layoutLegend(lg, 'south');
layoutTitle(t, sprintf(['RQ3. Sunlight-weighted Albedo of All Alaska (%s): observed %s, and temperature-scaled for ' ...
    'every year %s from the April-May warming of NEX-GDDP-CMIP6 (%+.3f per ^{\\circ}C)'], ss, by, fy, F.slope(j)), C);
saveFig(f, fullfile(E.figDir, 'fig19_temperature_sensitivity.pdf'));
end

function col = scenColour(scen)
% IPCC AR6 WGI colours of the SSP scenarios (fig19): one row per name, grey for any other name
k = {'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'};
c = [0 173 207; 23 60 102; 247 148 32; 231 29 37; 149 27 30] / 255;
col = repmat([0.5 0.5 0.5], numel(scen), 1);
for s = 1:numel(scen)
    q = find(strcmp(k, scen{s}));
    if ~isempty(q), col(s, :) = c(q, :); end
end
end

function drawAlbTempMap(AT, G, E)
% fig19b (E.figures 19; RQ3): a map of b, the change of the season's sunlight-weighted Albedo per degC of April-May air
% temperature (tab19), by ecoregion: every kept 1-km cell filled with the b of its ecoregion (ecoregion and land
% outlines as in the other maps), each ecoregion labelled with its name and b, with an asterisk where the E.albTempCI %
% range of b (resamples of the years) holds zero (a note at the bottom right says so, or that no range does); a text
% block with All Alaska's b, r and number of years and the E.whatIfYears changes of the scenarios with the smallest and
% the largest change (the mean of the models; tab19). Colours: the red half of C.div (red = darker, its grey
% midpoint = no change) from the 1-2-5 step at or below the most negative b to 0; all of C.div over +-the largest |b|
% (blue = brighter) if an ecoregion has b > 0. Grey (C.never): kept cells outside the nine ecoregions, and an ecoregion
% without a fit. Labels at the places of regionLabelPlaces (those of fig19c too): inside the ecoregion, or outside it
% with a leader line to its cell nearest the listed point.
if isempty(AT) || ~any(E.figures == 19), return; end
C = palette();  F = AT.fit;  nR = numel(G.order);  j0 = nR + 1;  nS = numel(AT.scen);
assert(isequal(reshape(AT.names, 1, []), [reshape(G.regNames(G.order), 1, []), {'All Alaska'}]), ...
    'fig19b: the regions of AT are not those of G');
assert(isfinite(F.slope(j0)), 'fig19b: no fit for All Alaska');
b = reshape(F.slope(1:nR), 1, []);  has = isfinite(b);
assert(any(has), 'fig19b: no ecoregion has a fit');
star = has & reshape(F.lo(1:nR), 1, []) <= 0 & reshape(F.hi(1:nR), 1, []) >= 0;   % the range holds zero
lim = niceLimits([b(has), 0]);
if lim(2) <= 0
    cl = [lim(1) 0];  cmap = cmapFrom(C.div(1:6, :));     % the red half of C.div: dark red ... grey midpoint
else
    cl = max(abs(lim)) * [-1 1];  cmap = C.divmap;
end
v = nan(numel(G.keep), 1);
for j = 1:nR, v(G.regC == G.order(j)) = b(j); end
f = newFig(7.4, 7.2, E);  t = tiles(f, 1, 1);  a = tileAt(t, 1, [1 1], C);
mapRGB(a, mapNoData(v, G, cmap, cl, C), G.xc, G.yc, G, C, false);
xl = get(a, 'XLim');  yl = get(a, 'YLim');
[P, Q] = regionLabelPlaces(G, xl, yl, 'fig19b');          % label places and leader-line ends (km); NaN: no line
for j = 1:nR
    c = G.order(j);  p = P(j, :);  q = Q(j, :);
    if all(isfinite(q)), plot(a, [p(1) q(1)], [p(2) q(2)], '-', 'Color', C.ink2, 'LineWidth', 0.6); end
    if has(j)
        val = sprintf('$%+.3f$', b(j));
        if star(j), val = [val '*']; end
    else
        val = 'no fit';
    end
    text(a, p(1), p(2), {latexSafe(G.short{c}), val}, 'FontSize', 8, 'Color', C.ink, 'HorizontalAlignment', ...
        'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'BackgroundColor', 'w', 'Margin', 1.5);
end
sl = cellfun(@scenLabel, AT.scen, 'UniformOutput', false);
dA = arrayfun(@(s) AT.S(s).sum.dAMean(j0), 1:nS);             % the mean of the models
assert(all(isfinite(dA)), 'fig19b: a scenario has no estimate for All Alaska');
[~, i1] = max(dA);  [~, i2] = min(dA);                       % the smallest and the largest darkening
mth = @(fmt, x) ['$', sprintf(fmt, x), '$'];  fy = AT.years;
blk = {'\textbf{Statewide}', sprintf('%s per $^\\circ$C ($r = %.2f$, %d years)', mth('%+.3f', F.slope(j0)), ...
    F.r(j0), F.n(j0)), sprintf('%d-%d mean, temperature-scaled:', fy(1), fy(2))};
if i1 == i2
    blk{end + 1} = sprintf('%s (%s)', mth('%+.3f', dA(i1)), sl{i1});
else
    blk{end + 1} = sprintf('%s (%s) to %s (%s)', mth('%+.3f', dA(i1)), sl{i1}, mth('%+.3f', dA(i2)), sl{i2});
end
text(a, xl(2) - 0.01 * diff(xl), yl(2) - 0.01 * diff(yl), blk, 'FontSize', 8.5, 'Color', C.ink, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Interpreter', 'latex', 'BackgroundColor', 'w', ...
    'EdgeColor', C.grid, 'Margin', 4);
ss = sprintf('%s - %s', doyLabel(AT.days(1)), doyLabel(AT.days(2)));
layoutColorbar(a, cmap, cl, 'Albedo change per $^\circ$C of spring air temperature', 'south', C);%; grey: no ecoregion
if any(star)
    nt = sprintf('*: its %g \\%% range (resamples of the years) holds 0', E.albTempCI);
else
    nt = sprintf('the %g \\%% range of every ecoregion (resamples of the years) excludes 0', E.albTempCI);
end
text(a, xl(2) - 0.01 * diff(xl), yl(1) + 0.01 * diff(yl), nt, 'FontSize', 8, 'Color', C.ink2, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Interpreter', 'latex');
layoutTitle(t, sprintf(['RQ3. Albedo change per ^{\\circ}C of April-May air temperature, by ecoregion (%d-%d; ' ...
    'sunlight-weighted Albedo of %s)'], AT.yb(1), AT.yb(end), ss), C);          % TeX interpreter: a degree sign
saveFig(f, fullfile(E.figDir, 'fig19b_sensitivity_map.pdf'));
end

function drawAlbTempChangeMaps(AT, G, E)
% fig19c (E.figures 19; RQ3; maps beside the table tab19c): the temperature-scaled Albedo of the season
% (sunlight-weighted) in E.whatIfYears by ecoregion, one map per scenario (in the order of E.whatIfScenarios, row by
% row; 2 x 2 for four) on one colour scale: every kept 1-km cell filled with the change dA of its ecoregion from the
% baseline mean (dAMean = b x the models' mean warming of the April-May blocks, the mean of the yearly means of fig19;
% tab19, tab19c), ecoregion and land outlines as in the other maps; each ecoregion labelled at the places of fig19b
% (regionLabelPlaces; the same leader lines) with its E.whatIfYears Albedo (aFutMean) over its change in brackets, with
% an asterisk where the E.albTempCI % range of every model with every resample of the years holds zero (the colour bar's
% label says so, or that no range does). Over each map the scenario, and All Alaska's Albedo, change, the models' mean
% warming of April-May and their number. Colours as in fig19b: the red half of C.div (red = darker, its grey
% midpoint = no change) from the 1-2-5 step at or below the most negative dA of all maps to 0; all of C.div over +-the
% largest |dA| (blue = brighter) if an ecoregion has dA > 0. Grey (C.never): kept cells outside the nine ecoregions, and
% an ecoregion without an estimate. The ranges themselves are in tab19 and tab19c
if isempty(AT) || ~any(E.figures == 19), return; end
C = palette();  nR = numel(G.order);  j0 = nR + 1;  nS = numel(AT.scen);  fy = AT.years;  yb = AT.yb;
assert(isequal(reshape(AT.names, 1, []), [reshape(G.regNames(G.order), 1, []), {'All Alaska'}]), ...
    'fig19c: the regions of AT are not those of G');
D = nan(nS, nR);  Af = nan(nS, nR);  star = false(nS, nR);  star0 = false(1, nS);   % dA, Albedo; range holds zero
for s = 1:nS
    m = AT.S(s).sum;
    assert(isfinite(m.dAMean(j0)), 'fig19c: %s has no estimate for All Alaska', AT.scen{s});
    D(s, :) = reshape(m.dAMean(1:nR), 1, []);  Af(s, :) = reshape(m.aFutMean(1:nR), 1, []);
    star(s, :) = isfinite(D(s, :)) & reshape(m.dALo(1:nR), 1, []) <= 0 & reshape(m.dAHi(1:nR), 1, []) >= 0;
    star0(s) = m.dALo(j0) <= 0 && m.dAHi(j0) >= 0;
end
has = isfinite(D);
assert(any(has(:)), 'fig19c: no ecoregion has an estimate');
assert(isequal(has, isfinite(Af)), 'fig19c: an ecoregion has a change without an Albedo or the reverse');
lim = niceLimits([reshape(D(has), 1, []), 0]);                  % one colour scale for every map
if lim(2) <= 0
    cl = [lim(1) 0];  cmap = cmapFrom(C.div(1:6, :));     % the red half of C.div: dark red ... grey midpoint
else
    cl = max(abs(lim)) * [-1 1];  cmap = C.divmap;
end
sl = cellfun(@scenLabel, AT.scen, 'UniformOutput', false);
mth = @(fmt, x) ['$', sprintf(fmt, x), '$'];  ast = {'', '*'};   % numbers in math mode: true minus signs
nc = ceil(sqrt(nS));  nr = ceil(nS / nc);                        % 2 x 2 for four scenarios
if nS == 1, f = newFig(7.4, 7.2, E); else, f = newFig(5 * nc, 4.6 * nr + 1.2, E); end   % one map: as fig19b
t = tiles(f, nr, nc);
for s = 1:nS
    a = tileAt(t, s, [1 1], C);
    v = nan(numel(G.keep), 1);
    for j = 1:nR, v(G.regC == G.order(j)) = D(s, j); end
    mapRGB(a, mapNoData(v, G, cmap, cl, C), G.xc, G.yc, G, C, false);
    if s == 1, [P, Q] = regionLabelPlaces(G, get(a, 'XLim'), get(a, 'YLim'), 'fig19c'); end   % as in fig19b
    for j = 1:nR
        if all(isfinite(Q(j, :)))
            plot(a, [P(j, 1) Q(j, 1)], [P(j, 2) Q(j, 2)], '-', 'Color', C.ink2, 'LineWidth', 0.6);
        end
        if has(s, j)
            val = {mth('%.3f', Af(s, j)), ['(', mth('%+.3f', D(s, j)), ast{star(s, j) + 1}, ')']};
        else
            val = 'no estimate';
        end
        text(a, P(j, 1), P(j, 2), val, 'FontSize', 8, 'Color', C.ink, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'BackgroundColor', 'w', 'Margin', 1.5);
    end
    m = AT.S(s).sum;
    panelTitle(a, sl{s}, C, sprintf('All Alaska: %s (%s%s) for %s of April-May warming (mean of %d models)', ...
        mth('%.3f', m.aFutMean(j0)), mth('%+.3f', m.dAMean(j0)), ast{star0(s) + 1}, ...
        sprintf('$%+.1f\\,^\\circ$C', m.dTMean(j0)), m.nModels(j0)));
end
if any(star(:)) || any(star0)
    nt = sprintf('*: its %g \\%% range (models $\\times$ resamples of the years) holds 0', E.albTempCI);
else
    nt = sprintf('the %g \\%% range of every change (models $\\times$ resamples of the years) excludes 0', ...
        E.albTempCI);
end
layoutColorbar(a, cmap, cl, {sprintf('Albedo change, %d-%d minus %d-%d (mean of the models); grey: no ecoregion', ...
    fy(1), fy(2), yb(1), yb(end)), sprintf('labels: the %d-%d Albedo and (its change)', fy(1), fy(2)), nt}, ...
    'south', C);
ss = sprintf('%s - %s', doyLabel(AT.days(1)), doyLabel(AT.days(2)));
layoutTitle(t, sprintf(['RQ3. Temperature-scaled sunlight-weighted Albedo (%s) of %d-%d and its change from %d-%d, ' ...
    'by ecoregion and scenario (NEX-GDDP-CMIP6)'], ss, fy(1), fy(2), yb(1), yb(end)), C);
saveFig(f, fullfile(E.figDir, 'fig19c_change_maps.pdf'));
end

function [P, Q] = regionLabelPlaces(G, xl, yl, tag)
% places of the ecoregion labels of the maps fig19b and fig19c (both label the same places), km of EPSG:3338 on a map
% with the limits xl, yl: P(j, :) the label of ecoregion G.order(j), Q(j, :) the end of its leader line (NaN: none). The
% places listed below were set on the 1-km grid of Prepared\grid_1km: inside the ecoregion, or outside it with a leader
% line to its cell nearest the listed point. Where a listed place does not fit the grid (outside the map, or more than 5
% cells from the ecoregion) the label sits at the ecoregion's cell nearest its centroid, and a note is printed (tag: the
% figure); a region not listed is labelled there too
nR = numel(G.order);  P = nan(nR, 2);  Q = nan(nR, 2);  dc = 5 * G.agg * G.cell / 1e3;         % 5 cells (km)
% label places (km, EPSG:3338): Level-2 name, label, point the leader line goes to ([]: the label sits inside)
LP = {'Arctic Tundra', [0 2110], [];  'Bering Tundra', [-470 1760], [];  'Bering Taiga', [-352 1290], []; ...
      'Intermontane Boreal', [235 1700], [];  'Alaska Range Transition', [163 1353], []; ...
      'Pacific Mountains Transition', [440 870], [446 1389];  'Coast Mountains Transition', [1180 1560], [560 1372]; ...
      'Coastal Rainforests', [1000 690], [1100 940];  'Aleutian Meadows', [-560 880], [-284 787]};
inMap = @(p) p(1) >= xl(1) && p(1) <= xl(2) && p(2) >= yl(1) && p(2) <= yl(2);
for j = 1:nR
    c = G.order(j);  k = find(strcmp(LP(:, 1), G.regNames{c}));
    p = [];  q = [];
    if ~isempty(k) && inMap(LP{k, 2})
        if isempty(LP{k, 3})
            [~, d] = regionCell(c, LP{k, 2}, G);
            if d <= dc, p = LP{k, 2}; end
        else
            [q, d] = regionCell(c, LP{k, 3}, G);
            if d <= dc, p = LP{k, 2}; else, q = []; end
        end
    end
    if isempty(p)
        p = regionCell(c, [], G);
        if ~isempty(k)
            fprintf('  %s: %s labelled at its cell nearest its centroid (its listed place does not fit)\n', tag, ...
                G.regNames{c});
        end
    end
    P(j, :) = p;
    if ~isempty(q), Q(j, :) = q; end
end
end

function [p, d] = regionCell(c, p0, G)
% centre (km) of the cell of ecoregion code c nearest the point p0 ([]: the ecoregion's centroid) and its distance (km)
[r, k] = find(G.ecoC == c);
assert(~isempty(r), 'ecoregion code %d has no cell', c);
xr = reshape(G.xc(k), [], 1);  yr = reshape(G.yc(r), [], 1);
if isempty(p0), p0 = [mean(xr), mean(yr)]; end
[d2, i] = min((xr - p0(1)) .^ 2 + (yr - p0(2)) .^ 2);
p = [xr(i), yr(i)];  d = sqrt(d2);
end

function writeAlbTempTables(AT, G, E)
% tab19 (RQ3): one row per region and scenario: the season's sunlight-weighted Albedo and April-May temperature, the
% sensitivity (fit, intervals, checks) and the temperature-scaled estimate for E.whatIfYears (the median of the
% models; the period and the mean of the models in the last columns); tab19b: every region and year (the series of
% fig19, the retrieval mix); tab19c: the table for the paper (below); tab19d and tab19e: every year of E.whatIfYears
% (below). Every run with E.albTemp
if isempty(AT), return; end
full = [G.regNames(G.order), {'All Alaska'}];  F = AT.fit;  B = AT.blk;  ci = E.albTempCI;
cl = @(s) {sprintf('%s_ci%g_low%s', s, ci, ''), sprintf('%s_ci%g_high%s', s, ci, '')};
hdr = [{'region', 'scenario', 'years', 'season_first_doy', 'season_last_doy', 'sunlight_share_of_year', 't_first_doy', ...
    't_last_doy', 'albedo_eff_mean', 'albedo_eff_sd', 't_mean_degC', 't_sd_K', 'slope_per_K', 'slope_se'}, ...
    cl('slope'), {'r', 'p', 'slope_detrended_per_K', 'lag1_autocorr_residuals', 'blockwise_slope_per_K'}, ...
    cl('blockwise_slope'), {'models', 'warming_K'}, cl('warming'), {'warming_in_sd', 'baseline_years_as_warm', ...
    't_future_degC', 'albedo_change'}, cl('albedo_change'), cl('albedo_change_models'), ...
    cl('albedo_change_resamples'), {'albedo_eff_future'}, cl('albedo_eff_future'), {'extra_absorbed_MJ_m2_yr'}, ...
    cl('extra_absorbed_MJ_m2_yr'), {'extra_absorbed_W_m2'}, cl('extra_absorbed_W_m2'), ...
    {'blockwise_albedo_change'}, cl('blockwise_albedo_change'), cl('slope_t'), {'slope_leave_one_out_min', ...
    'slope_leave_one_out_max', 'year_trend_per_decade', 'year_trend_se', 'year_trend_p', ...
    'slope_with_magnitude_share_per_K', 'slope_with_magnitude_share_se', 'slope_with_valid_share_per_K', ...
    'slope_with_valid_share_se', 'slope_part_before_t_window', 'slope_part_t_window', 'slope_part_after_t_window', ...
    'future_first_year', 'future_last_year', 'warming_mean_K', 'albedo_change_mean', 'albedo_eff_future_mean', ...
    'extra_absorbed_W_m2_mean'}];
rows = cell(10 * numel(AT.S), numel(hdr));  q = 0;
for j = 1:10
    for s = 1:numel(AT.S)
        m = AT.S(s).sum;  q = q + 1;
        rows(q, :) = {full{j}, AT.scen{s}, F.n(j), AT.days(1), AT.days(2), AT.share(j), AT.tDays(1), AT.tDays(2), ...
            F.meanA(j), F.sdA(j), F.meanT(j), F.sdT(j), F.slope(j), F.se(j), F.lo(j), F.hi(j), F.r(j), F.p(j), ...
            F.slopeDet(j), F.r1(j), B.U(j), B.ULo(j), B.UHi(j), m.nModels(j), m.dT(j), m.dTLo(j), m.dTHi(j), ...
            m.dTsd(j), m.asWarm(j), m.tFut(j), m.dA(j), m.dALo(j), m.dAHi(j), m.dAmLo(j), m.dAmHi(j), m.dAbLo(j), ...
            m.dAbHi(j), m.aFut(j), m.aFutLo(j), m.aFutHi(j), m.dQ(j), m.dQLo(j), m.dQHi(j), m.dQW(j), m.dQWLo(j), ...
            m.dQWHi(j), m.dAblk(j), m.dAblkLo(j), m.dAblkHi(j), F.tLo(j), F.tHi(j), F.looMin(j), F.looMax(j), ...
            F.yearTrend(j), F.yearSE(j), F.yearP(j), F.slopeShare(j), F.slopeShareSE(j), F.slopeValid(j), ...
            F.slopeValidSE(j), F.part(j, 1), F.part(j, 2), F.part(j, 3), AT.years(1), AT.years(2), m.dTMean(j), ...
            m.dAMean(j), m.aFutMean(j), m.dQWMean(j)};
    end
end
writeCsv(fullfile(E.tabDir, 'tab19_albedo_temperature.csv'), hdr, rows);
nb = numel(AT.yb);  rows = cell(10 * nb, 9);  q = 0;
for j = 1:10
    for i = 1:nb
        q = q + 1;
        fit = F.icp(j) + F.slope(j) * AT.t(j, i);
        rows(q, :) = {full{j}, AT.yb(i), AT.alb(j, i), AT.t(j, i), double(F.used(j, i)), fit, AT.alb(j, i) - fit, ...
            AT.qaShare(j, i), AT.qaValid(j, i)};
    end
end
writeCsv(fullfile(E.tabDir, 'tab19b_albedo_temperature_years.csv'), {'region', 'year', 'albedo_eff', 't_mean_degC', ...
    'in_fit', 'albedo_eff_from_t', 'residual', 'magnitude_share', 'valid_share'}, rows);
% tab19c (a table for the paper): one row per region, each cell as text 'value (low to high)' with E.albTempCI % ranges:
% b per degC (its range from the resamples of the years) and each scenario's temperature-scaled change of the
% E.whatIfYears mean Albedo (the mean of the models; the range of every model with every resample); the numbers of
% tab19, rounded; the header gives each scenario's period and models (a range if the regions differ)
nS = numel(AT.S);  hdr = cell(1, 2 + nS);  rows = cell(10, 2 + nS);
hdr(1:2) = {'region', sprintf('Albedo change per degC %d-%d', AT.yb(1), AT.yb(end))};
for s = 1:nS
    m = AT.S(s).sum;  nm = m.nModels(isfinite(m.dAMean));
    if isempty(nm), mt = 'no models'; elseif min(nm) == max(nm), mt = sprintf('%d models', nm(1));
    else, mt = sprintf('%d-%d models', min(nm), max(nm)); end
    hdr{2 + s} = sprintf('%s %d-%d (%s)', scenLabel(AT.scen{s}), AT.years(1), AT.years(2), mt);
end
medRange = @(fmt, v) sprintf([fmt ' (' fmt ' to ' fmt ')'], v);          % v = [value low high]
for j = 1:10
    rows(j, 1:2) = {full{j}, ''};
    if isfinite(F.slope(j)), rows{j, 2} = medRange('%+.4f', [F.slope(j), F.lo(j), F.hi(j)]); end
    for s = 1:nS
        m = AT.S(s).sum;  rows{j, 2 + s} = '';
        if isfinite(m.dAMean(j)), rows{j, 2 + s} = medRange('%+.3f', [m.dAMean(j), m.dALo(j), m.dAHi(j)]); end
    end
end
writeCsv(fullfile(E.tabDir, 'tab19c_albedo_temperature_paper.csv'), hdr, rows);
% tab19d: one row per region, scenario and year of E.whatIfYears: the models' warming of April-May (mean,
% median, E.albTempCI % range) and the temperature-scaled Albedo A(y) = the baseline mean + b dT(y) (mean, median,
% quartiles, E.albTempCI % range of the models; the boxes of fig19), its change (mean) and how many models are warmer
% than the warmest fitted year; tab19e: every model's warming and A(y). Empty: no fit
nYf = numel(AT.yf);
hdr = [{'region', 'scenario', 'year', 'models', 'warming_mean_K', 'warming_median_K'}, cl('warming'), ...
    {'albedo_eff_mean', 'albedo_eff_median', 'albedo_eff_q25', 'albedo_eff_q75'}, cl('albedo_eff'), ...
    {'albedo_change_mean', 'models_warmer_than_baseline_years'}];
rows = cell(10 * nS * nYf, numel(hdr));  q = 0;
for j = 1:10
    for s = 1:nS
        y = AT.S(s).yr;  nm = numel(AT.S(s).models);
        for i = 1:nYf
            q = q + 1;
            rows(q, :) = {full{j}, AT.scen{s}, AT.yf(i), nm, y.dTMean(j, i), y.dTMed(j, i), y.dTLo(j, i), ...
                y.dTHi(j, i), y.aMean(j, i), y.aMed(j, i), y.aQ1(j, i), y.aQ3(j, i), y.aLo(j, i), y.aHi(j, i), ...
                y.dAMean(j, i), y.nWarmer(j, i)};
        end
    end
end
writeCsv(fullfile(E.tabDir, 'tab19d_albedo_yearly.csv'), hdr, rows);
nM = arrayfun(@(s) numel(AT.S(s).models), 1:nS);
rows = cell(10 * sum(nM) * nYf, 6);  q = 0;
for j = 1:10
    for s = 1:nS
        S = AT.S(s);
        for k = 1:nM(s)
            for i = 1:nYf
                q = q + 1;
                rows(q, :) = {full{j}, AT.scen{s}, S.models{k}, AT.yf(i), S.dTy(k, j, i), S.aY(k, j, i)};
            end
        end
    end
end
writeCsv(fullfile(E.tabDir, 'tab19e_albedo_yearly_models.csv'), {'region', 'scenario', 'model', 'year', 'warming_K', ...
    'albedo_eff'}, rows);
fprintf('  wrote tab19, tab19b, tab19c, tab19d and tab19e to %s\n', E.tabDir);
end

% ----------------------------- statistics (no toolbox) -----------------------------

function beta = lsq(Xd, y, what)
% least squares through a QR decomposition; error if the design is (numerically) rank deficient (what: the fit, named
% in the error message)
[Q, Rq] = qr(Xd, 0);  d = abs(diag(Rq));
assert(min(d) > 1e-9 * max(d), '%s: rank-deficient design (%d x %d): a covariate does not vary over the years', what, ...
    size(Xd, 1), size(Xd, 2));
beta = Rq \ (Q' * y);
end

function t = tQuantile(ci, df)
% two-sided ci-% critical value of Student's t with df degrees of freedom (its (1 + ci/100)/2 quantile), from the
% incomplete beta function: P(|T| > t) = I(df / (df + t^2); df/2, 1/2)
x = betaincinv(1 - ci / 100, df / 2, 0.5);
t = sqrt(df * (1 - x) / x);
end

function B = to365(A, y)
% columns = actual day of year 1..366 -> 365-day calendar (29 Feb dropped)
if mod(y, 4) == 0 && (mod(y, 100) ~= 0 || mod(y, 400) == 0)
    B = A(:, [1:59, 61:366]);
else
    B = A(:, 1:365);
end
end

function m = nanMean(X, dim)
% mean over dimension dim ignoring NaN (NaN where none)
ok = isfinite(X);  X(~ok) = 0;
n = sum(ok, dim);
m = sum(X, dim) ./ n;  m(n == 0) = NaN;
end

function s = nanStd(X, dim)
% sample SD (n - 1) over dimension dim ignoring NaN (NaN with < 2 values)
ok = isfinite(X);  n = sum(ok, dim);
mu = nanMean(X, dim);
d = X - mu;  d(~ok) = 0;
s = sqrt(sum(d.^2, dim) ./ (n - 1));  s(n < 2) = NaN;
end

function out = nanMovMean(x, half, minN)
% centred moving mean over +-half samples, NaN-aware; NaN with < minN values
x = x(:)';  ok = isfinite(x);  v = x;  v(~ok) = 0;
k = ones(1, 2 * half + 1);
s = conv(v, k, 'same');  n = conv(double(ok), k, 'same');
out = s ./ max(n, 1);  out(n < minN) = NaN;
end

function swIn = shortwave(W, yb, E, order)
% Daily climatology of the incoming shortwave (ERA5-Land SSRD, MJ m-2 d-1; the weights of stage 2h). W: the sums of
% passShortwave [11 x 366 x years]; yb: the years. Rows: the 9 regions (order), then statewide; 365-day calendar.
% A region's value of a day = the pooled mean over its ERA5-Land pixels; swIn.mean, .sd [10 x 365]: the mean and SD
% over the years (>= E.minCycleYears), .meanS, .sdS: 7-day smoothed (climOf)
k2 = 0.002;                                                % SSRD raw -> MJ m-2 d-1
nb = numel(yb);
rowsOf = @(X) [X(order + 1, :); sum(X(2:11, :), 1)];
in = nan(10, 365, nb);
for i = 1:nb
    nA = to365(rowsOf(W.nAll(:, :, i)), yb(i));  iA = to365(rowsOf(W.inAll(:, :, i)), yb(i));
    v = iA ./ max(nA, 1) * k2;  v(nA == 0) = NaN;  in(:, :, i) = v;
end
swIn = climOf(in, E);
end

function C = climOf(X, E)
% [10 x 365 x years] -> daily mean and SD over the years (>= E.minCycleYears values), raw and 7-day smoothed
nY = size(X, 3);  C = struct('mean', nan(10, 365), 'sd', nan(10, 365), 'meanS', nan(10, 365), 'sdS', nan(10, 365));
for j = 1:10
    Y = reshape(X(j, :, :), 365, nY);
    n = sum(isfinite(Y), 2);  mu = nanMean(Y, 2);  sd = nanStd(Y, 2);
    mu(n < E.minCycleYears) = NaN;  sd(n < E.minCycleYears) = NaN;
    C.mean(j, :) = mu;  C.sd(j, :) = sd;
    C.meanS(j, :) = nanMovMean(mu, 3, 4);  C.sdS(j, :) = nanMovMean(sd, 3, 4);
end
end

function clim = blockNormal(M, minYears)
% cells x blocks x years (8-day means) -> cells x blocks: the cell's mean of each block over
% the years (summed in year order, as the Python twin; NaN with < minYears values)
[nK, nc, nY] = size(M);
s = zeros(nK, nc);  n = zeros(nK, nc);
for y = 1:nY
    x = double(M(:, :, y));  f = isfinite(x);  x(~f) = 0;
    s = s + x;  n = n + f;
end
clim = s ./ max(n, 1);  clim(n < minYears) = NaN;
end

function Z = blockAnomalies(M, minYears)
% cells x blocks x years -> cells x (blocks * years) anomalies (value - blockNormal),
% block index fastest
M = double(M);  [nK, nc, nY] = size(M);
Z = reshape(M - blockNormal(M, minYears), nK, nc * nY);
end

function [o, lin, first, last] = rowRuns(X)
% sort every row (NaN last; sort is stable) and find, for each sorted position, the first
% and last position of its run of equal values (NaN ~= NaN: a run of its own)
[nr, nc] = size(X);
[s, o] = sort(X, 2);
lin = sub2ind([nr nc], repmat((1:nr)', 1, nc), o);                 % original element of each sorted slot
pos = repmat(1:nc, nr, 1);
st = [true(nr, 1), s(:, 2:end) ~= s(:, 1:end-1)];
en = [s(:, 2:end) ~= s(:, 1:end-1), true(nr, 1)];
first = cummax(st .* pos, 2);
last = nc + 1 - fliplr(cummax(fliplr(en .* (nc + 1 - pos)), 2));
end

function R = rowRankAvg(X)
% average ranks (1..n) of the finite values of every row; NaN stays NaN
[~, lin, first, last] = rowRuns(X);
R = zeros(size(X));  R(lin) = (first + last) / 2;
R(~isfinite(X)) = NaN;
end

function [r, n] = rowSpearman(A, B, minN)
% Spearman rho of every row pair over the columns where both are finite (average ranks);
% NaN with < minN pairs or no variance. Also the number of pairs.
m = isfinite(A) & isfinite(B);  A(~m) = NaN;  B(~m) = NaN;
ra = rowRankAvg(A);  rb = rowRankAvg(B);
n = sum(m, 2);  mu = (n + 1) / 2;
da = ra - mu;  db = rb - mu;  da(~m) = 0;  db(~m) = 0;
r = sum(da .* db, 2) ./ sqrt(sum(da .^ 2, 2) .* sum(db .^ 2, 2));
r(n < minN) = NaN;
end

function xi = rowXi(X, Y, K, minN)
% Chatterjee's xi(x -> y) of every row over the columns where X and Y are finite, as
% xiCorr: order by x, then key K, then y; r = #{y_j <= y_i}, l = #{y_j >= y_i}
m = isfinite(X) & isfinite(Y);  [nr, nc] = size(X);
X(~m) = NaN;  Y(~m) = NaN;  n = sum(m, 2);
[~, lin, first, last] = rowRuns(Y);
rr = zeros(nr, nc);  ll = zeros(nr, nc);
rr(lin) = last;  ll(lin) = n - first + 1;
ri = repmat((1:nr)', 1, nc);
[~, o] = sort(Y, 2);  lo = sub2ind([nr nc], ri, o);                % lexicographic (x, key, y):
[~, o2] = sort(K(lo), 2);  lo = lo(sub2ind([nr nc], ri, o2));      % stable sorts by y, key, x
[~, o3] = sort(X(lo), 2);  lo = lo(sub2ind([nr nc], ri, o3));
rs = rr(lo);  ms = m(lo);                                          % valid ones first
num = sum(abs(diff(rs, 1, 2)) .* (ms(:, 2:end) & ms(:, 1:end-1)), 2);
den = sum(ll .* (n - ll) .* m, 2);
xi = 1 - n .* num ./ (2 * den);
xi(n < minN) = NaN;
end

function n = northness(slope, aspect)
% sin(slope) * cos(aspect) (degrees; aspect clockwise from north): +1 steep north-facing, -1 steep
% south-facing, 0 flat or east/west-facing. Rounded to 1e-4 (half up) so that MATLAB and the Python
% twin, whose sin/cos may differ in the last bit, give identical values and ranks.
n = floor(sin(slope * (pi / 180)) .* cos(aspect * (pi / 180)) * 1e4 + 0.5) / 1e4;
end

function k = tieKey(pix, year, doy)
% fixed pseudo-random tie-break key of a pixel-day, identical in the Python twin: a
% linear seed of (row-major 1-km index, year, day of year) mod 2^32 mixed by 4 rounds
% of Marsaglia's xorshift32 (shifts 13, 17, 5; overflowing bits dropped)
s = uint32(mod(double(pix) * 7919 + double(year) * 104729 + double(doy) * 1299709, 2^32));
for it = 1:4
    s = bitxor(s, bitshift(s, 13));  s = bitxor(s, bitshift(s, -17));  s = bitxor(s, bitshift(s, 5));
end
k = double(s);
end

function xi = xiCorr(x, y, key)
% Chatterjee's xi(x -> y) (Chatterjee 2021, JASA; formula with ties in y): how far y is a
% function of x, 0 (independent) .. 1. Rows ordered by x, ties in x broken by key, then y.
n = numel(x);
[~, o] = sortrows([x(:), key(:), y(:)]);
[~, ~, ic] = unique(y(:));
cntU = accumarray(ic, 1);  cum = cumsum(cntU);
rAll = cum(ic);                                        % #{j: y_j <= y_i}
lAll = n - (cum(ic) - cntU(ic));                       % #{j: y_j >= y_i}
r = rAll(o);  l = lAll(o);
xi = 1 - n * sum(abs(diff(r))) / (2 * sum(l .* (n - l)));
end

function r = rankavg(x)
% ranks 1..n, ties get their average rank
x = x(:);  n = numel(x);
[xs, ord] = sort(x);
st = [true; diff(xs) ~= 0];
first = find(st);  last = [first(2:end) - 1; n];
avg = (first + last) / 2;
r = zeros(n, 1);  r(ord) = avg(cumsum(st));
end

function r = pearson(x, y)
x = x(:) - mean(x(:));  y = y(:) - mean(y(:));
r = sum(x .* y) / sqrt(sum(x .* x) * sum(y .* y));
end

function [Rm, N] = spearmanMatrix(V)
% pairwise-complete Spearman matrix of the columns of V and the pair counts
k = size(V, 2);  Rm = eye(k);  N = zeros(k);  fin = isfinite(V);
for a = 1:k
    N(a, a) = nnz(fin(:, a));
    for b = a + 1:k
        m = fin(:, a) & fin(:, b);
        N(a, b) = nnz(m);  N(b, a) = N(a, b);
        if N(a, b) > 2
            Rm(a, b) = pearson(rankavg(V(m, a)), rankavg(V(m, b)));
        else
            Rm(a, b) = NaN;
        end
        Rm(b, a) = Rm(a, b);
    end
end
end

function q = pctl(x, p)
% percentile, linear interpolation between order statistics (numpy 'linear'), NaN ignored
x = sort(x(isfinite(x)));  n = numel(x);
if n == 0, q = NaN; return; end
pos = p / 100 * (n - 1);  lo = floor(pos);  hi = min(lo + 1, n - 1);
q = x(lo + 1) + (pos - lo) * (x(hi + 1) - x(lo + 1));
end

function c = monthNames(m)
nm = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
c = nm(m);
end

function printSummary(R, E)
fprintf('Never-observed AOI pixels (%s-%s, baseline): %.1f %%\n', R.mname{1}, R.mname{end}, 100 * R.neverShare);
for si = 1:numel(R.blkSeasons)
    c = {};
    for di = 1:numel(R.blkDrivers)
        r = squeeze(R.blkRho(si, di, :));  r = r(isfinite(r));  mr = NaN;
        if ~isempty(r), mr = median(r); end
        c = [c, {R.blkDrivers{di}, mr, numel(r)}];                         %#ok<AGROW>
    end
    fprintf('%s driver anomalies, median rho (cells): %s\n', R.blkSeasons{si}, sprintf('%s %+.2f (%d)  ', c{:}));
end
c = [arrayfun(@(q) sprintf('%d-%d', R.periods(q, 1), R.periods(q, 2)), 1:size(R.periods, 1), 'UniformOutput', false); ...
     num2cell(R.perReg(10, :) - R.perRegBase(10))];
fprintf('Five-year Feb-Oct mean Albedo minus the baseline, statewide (%d cells with every period): %s\n', R.perRegN(10), ...
    sprintf('%s %+.3f  ', c{:}));
fprintf('Cells per permafrost class (none, isolated, sporadic, discontinuous, continuous): %s\n', mat2str(R.pfCells));
for s = 1:3
    c = [R.vars(2:end); num2cell(R.rho{s}(2:end)); num2cell(R.xi{s}(2:end))];
    fprintf('%s rho / xi: %s\n', E.seasons{s, 1}, sprintf('%s %+.2f / %.2f  ', c{:}));
end
end

% ----------------------------- figures -----------------------------

function drawFigures(R, G, E, clr_yrs, X)
C = palette();
yrs = R.years;  ib = R.ib;  yb = yrs(ib);  nb = numel(yb);
bs = sprintf('%d-%d', R.base(1), R.base(2));
names = [G.short(G.order), {'All Alaska'}];
P = [10 1:9];                                   % panel order: statewide first, then the ecoregions
seasN = E.seasons(:, 1)';
xd = 1:365;  win = xd >= 32 & xd <= 304;
monTick = {'XTick', [32 91 152 213 274], 'XTickLabel', {'Feb', 'Apr', 'Jun', 'Aug', 'Oct'}};

% ---------- fig01: mean annual cycle with every baseline year ----------
if any(E.figures == 1)
    cyc = struct('yb', yb, 'names', {names}, 'P', P, 'xd', xd, 'win', win, 'monTick', {monTick}, 'bs', bs);
    if isempty(E.yearGroups)                                 % one colour per year: light (first) to dark (last)
        ycm = cmapFrom(clr_yrs);
        cyc.ycol = ycm(round(linspace(1, size(ycm, 1), nb)), :);  cyc.grp = [];
    else                                                     % one colour per group of years (E.yearGroups)
        cyc.grp = yearGroups(E.yearGroups, yb, clr_yrs);
    end
    % annual metrics of the legends: Albedo from X.met (tab11)
    met = struct('val', X.met, 'name', {X.metName}, 'unit', {X.metUnit}, 'vunit', '', 'fmt', '%.3f', ...
        'show', groupMetrics(E, 'Albedo'));
    drawCycleFig(R.cycYear, R.cycMean, met, [0 0.8], 'Blue Sky Albedo', 'Albedo', 'regional daily mean, 7-day smoothed', ...
        fullfile(E.figDir, 'fig01_mean_annual_cycle.pdf'), cyc, C, E);

    % ---------- fig01 of every driver: the same annual cycle of NDSI, NDVI, t2m, SWE and snowfall ----------
    assert(isequal(X.drv, {'ndsi', 'ndvi', 't2m', 'swe', 'snowfall'}), 'fig01: driver labels for %s', strjoin(X.drv, ', '));
    dlab = {'NDSI snow cover (0-1)', 'NDVI', '2-m air temperature ($^\circ$C)', 'SWE (mm)', 'snowfall (mm d$^{-1}$)'};
    dunit = {'', '', 'deg C', 'mm', 'mm/d'};  dfmt = {'%.3f', '%.3f', '%.1f', '%.0f', '%.2f'};   % legend metrics
    for v = 1:numel(X.drv)
        Y = X.drvCycYear(:, win, :, v);  yl = niceLimits(double(Y(isfinite(Y))));
        met = struct('val', [], 'name', {X.metName}, 'unit', {X.metUnit}, 'vunit', dunit{v}, 'fmt', dfmt{v}, ...
            'show', groupMetrics(E, X.drv{v}));
        if ~isempty(E.yearGroups)                            % the annual metrics of the driver's own curves, as
            met.val = nan(10, nb, numel(X.metName));         % those of Albedo (same windows and rules)
            for j = 1:10
                for i = 1:nb
                    met.val(j, i, :) = cycleMetrics(double(X.drvCycYear(j, :, i, v)), X.drvCycMean(j, :, v), E.metricMinFrac);
                end
            end
        end
        drawCycleFig(X.drvCycYear(:, :, :, v), X.drvCycMean(:, :, v), met, yl, dlab{v}, X.drvName{v}, ...
            'regional daily mean over the valid 1-km pixels, 7-day smoothed', ...
            fullfile(E.figDir, sprintf('fig01_driver_%s_annual_cycle.pdf', X.drv{v})), cyc, C, E);
    end
end

% ---------- fig16: Albedo and NDSI snow cover annual cycles in the same panels ----------
if any(E.figures == 16)
    assert(~isempty(E.yearGroups), 'fig16 draws the groups of years of E.yearGroups: set E.yearGroups');
    assert(isequal(size(E.bothYLim), [1 2]) && diff(E.bothYLim) > 0, ['E.bothYLim: [lo hi] of the axis shared by ' ...
        'Albedo and NDSI snow cover']);
    cyc = struct('yb', yb, 'P', P, 'xd', xd, 'win', win, 'bs', bs, 'grp', yearGroups(E.yearGroups, yb, clr_yrs));
    cyc.names = names;  cyc.monTick = monTick;
    v = find(strcmp(X.drv, 'ndsi'));
    assert(isscalar(v) && size(X.met, 1) == 10 && size(X.met, 2) == nb, 'fig16: NDSI curves and Albedo metrics expected');
    hmA = X.met(:, :, 7);               % half-melt date of every region and year: Albedo (as the fig01 legends)
    hmN = nan(10, nb);                  % ... and NDSI snow cover (as the legends of fig01_driver_ndsi_annual_cycle)
    for j = 1:10
        for i = 1:nb
            m = cycleMetrics(double(X.drvCycYear(j, :, i, v)), X.drvCycMean(j, :, v), E.metricMinFrac);
            hmN(j, i) = m(7);
        end
    end
    drawCycleBoth(R.cycYear, R.cycMean, X.drvCycYear(:, :, :, v), X.drvCycMean(:, :, v), hmA, hmN, ...
        fullfile(E.figDir, 'fig16_Albedo_ndsi_annual_cycle.pdf'), cyc, C, E);
end

% ---------- fig17: whole-period daily series of Albedo, NDSI, SWE and air temperature ----------
if any(E.figures == 17)
    iv = zeros(1, 3);  dn = {'ndsi', 't2m', 'swe'};
    for k = 1:3
        m = find(strcmp(X.drv, dn{k}));
        assert(isscalar(m), 'fig17: driver %s not found in X.drv', dn{k});
        iv(k) = m;
    end
    assert(size(R.cycYear, 3) == nb && size(X.drvCycYear, 3) == nb, 'fig17: curves of the %d baseline years expected', nb);
    assert(numel(E.seriesTLim) == 2 && diff(E.seriesTLim) > 0, 'E.seriesTLim: [lo hi] in degC');
    drawSeriesFig(R.cycYear, X.drvCycYear(:, :, :, iv), yb, names, P, ...
        fullfile(E.figDir, 'fig17_Albedo_drivers_timeseries.pdf'), C, E);
end

% ---------- fig02: maps of Albedo vs driver anomalies (Spearman rho per cell) ----------
if any(E.figures == 2)
    dl = {'t2m', 'SWE', 'snowfall', 'NDSI', 'NDVI'};
    sl = containers.Map({'All', 'MAM', 'JJA', 'SON'}, {'All year', 'MAM', 'JJA', 'SON'});
    nS = numel(R.blkSeasons);  nD = numel(R.blkDrivers);
    f = newFig(18, 4 + 3.1 * nS, E);
    t = tiles(f, nS, nD);
    for si = 1:nS
        sn = R.blkSeasons{si};  if isKey(sl, sn), sn = sl(sn); end
        for di = 1:nD
            r = squeeze(R.blkRho(si, di, :));
            rG = toGridC(r, G);
            rgb = toRGB(rG, C.divmap, [-0.8 0.8], [1 1 1]);
            gap = false(G.HC, G.WC);  gap(G.keep) = ~isfinite(r(:));       % cells with too few pairs
            % Replace each RGB channel at missing cells with the "never observed" color.
            for ch = 1:3, L = rgb(:, :, ch);  L(gap) = C.never(ch);  rgb(:, :, ch) = L; end
            a = tileAx(t, (si - 1) * nD + di, [1 1], C);
            mapRGB(a, rgb, G.xc, G.yc, G, C, false);
            rr = r(isfinite(r));
            if isempty(rr), sub = 'too few pairs'; else, sub = sprintf('median: %+.2f', median(rr)); end
            panelTitle(a, sprintf('%s: %s', sn, dl{di}), C, sub);
        end
    end
    layoutColorbar(a, C.divmap, [-0.8 0.8], ...
        'Spearman correlation coefficient $\rho$', 'south', C);
    layoutTitle(t, sprintf(['Where Albedo anomalies follow which driver: %d-km cells, 8-day blocks, %s ' ...
        '(All year: every block; grey: < %d pairs)'], G.cellKm, bs, E.blkMinPairs), C);
    saveFig(f, fullfile(E.figDir, 'fig02_driver_anomaly_maps.pdf'));
end

% ---------- fig03: Spearman rho and Chatterjee xi of Albedo vs each covariate, by region ----------
% rows: covariates; columns: All Alaska, then the 9 ecoregions (as fig01); one heat map per season and measure
if any(E.figures == 3)
    cv = R.vars(2:end);  nv = numel(cv);  nDay = 5;                % the first nDay covariates: same pixel-day
    dp = X.dep;  colsR = P;  nR = numel(colsR);                    % dep rows: ecoregions (G.order), then statewide
    f = newFig(22, 4.2 + 0.55 * 2 * nv, E);
    t = tiles(f, 2, 3);
    spec = {dp.rho, C.divmap, [-1 1], 'Spearman \rho', [-1 -0.5 0 0.5 1]; ...
            dp.xi, C.seq, [0 1], 'Chatterjee \xi(covariate \rightarrow Albedo)', 0:0.25:1};
    for q = 1:2
        for s3 = 1:3
            a = tileAt(t, (q - 1) * 3 + s3, [1 1], C);
            M = spec{q, 1}(colsR, :, s3)';                         % covariates x regions
            image(a, 1:nR, 1:nv, toRGB(M, spec{q, 2}, spec{q, 3}, [1 1 1]));
            set(a, 'YDir', 'reverse', 'XLim', [0.5 nR + 0.5], 'YLim', [0.5 nv + 0.5], 'XTick', 1:nR, 'YTick', 1:nv, ...
                'TickLength', [0 0]);
            if s3 == 1, set(a, 'YTickLabel', cv); else, set(a, 'YTickLabel', {}); end
            if q == 2, set(a, 'XTickLabel', names(colsR), 'XTickLabelRotation', 40); else, set(a, 'XTickLabel', {}); end
            for i = 1:nv
                for k = 1:nR
                    v = M(i, k);
                    if ~isfinite(v), continue; end
                    if q == 1, str = sprintf('%+.2f', v); dk = abs(v) > 0.6; else, str = sprintf('%.2f', v); dk = v > 0.55; end
                    text(a, k, i, str, 'HorizontalAlignment', 'center', 'FontSize', 7.5, 'Color', textOn(dk, C));
                end
            end
            plot(a, [0.5 nR + 0.5], (nDay + 0.5) * [1 1], '-', 'Color', C.ink, 'LineWidth', 1);
            plot(a, 1.5 * [1 1], [0.5 nv + 0.5], '-', 'Color', C.ink, 'LineWidth', 1);   % statewide | ecoregions
            panelTitle(a, sprintf('%s  %s', seasN{s3}, spec{q, 4}), C);
            if s3 == 3, layoutColorbar(a, spec{q, 2}, spec{q, 3}, '', 'axes', C, spec{q, 5}); end
        end
    end
    layoutTitle(t, sprintf(['Albedo against each covariate within each ecoregion (sample pixels of the region, %s; above the ' ...
        'line: same pixel-day; below: pixel mean Albedo of the season; permafrost: extent class 0 none .. 4 continuous; ' ...
        'blank: < %d pixel-days or pixels, or no variation)'], bs, E.regMinN), C);
    saveFig(f, fullfile(E.figDir, 'fig03_dependence_rho_xi.pdf'));
end

% ---------- fig07: five-year mean Albedo (Feb-Oct) and its difference from the baseline normal ----------
if any(E.figures == 7)
    Pp = R.periods;  nP = size(Pp, 1);
    lastDay = datenum(R.yearTable(end).last, 'yyyy-mm-dd');            % last day read (partial last year?)
    f = newFig(3.1 * nP + 1.4, 7.8, E);
    t = tiles(f, 2, 1);
    t1 = tilesIn(t, 1, nP, 1, [1 1]);  t2 = tilesIn(t, 1, nP, 2, [1 1]);
    for q = 1:nP
        lab = sprintf('%d-%d', Pp(q, 1), Pp(q, 2));
        note = {};
        if Pp(q, 2) >= max(yrs) && lastDay < datenum(max(yrs), 12, 31)
            note{end + 1} = sprintf('to %s', R.yearTable(end).last);                                        %#ok<AGROW>
        end
        if R.perYears(q) < Pp(q, 2) - Pp(q, 1) + 1, note{end + 1} = sprintf('%d years read', R.perYears(q)); end %#ok<AGROW>
        if ~isempty(note), lab = sprintf('%s (%s)', lab, strjoin(note, '; ')); end
        a = tileAx(t1, q, [1 1], C);
        mapRGB(a, mapNoData(double(R.perMean(:, q)), G, C.seq, [0.1 0.7], C), G.xc, G.yc, G, C, false);
        panelTitle(a, lab, C, sprintf('statewide: %.3f', R.perReg(10, q)));
        a2 = tileAx(t2, q, [1 1], C);
        mapRGB(a2, mapNoData(double(R.perMean(:, q) - R.perBase), G, C.divmap, [-0.1 0.1], C), G.xc, G.yc, G, C, false);
        panelTitle(a2, lab, C, sprintf('statewide: %+.3f', R.perReg(10, q) - R.perRegBase(10)));
    end
    layoutColorbar(a, C.seq, [0.1 0.7], 'Feb-Oct mean Albedo', 'east', C);
    layoutColorbar(a2, C.divmap, [-0.1 0.1], 'Anomaly', 'east', C);
    layoutTitle(t, sprintf(['Five-year mean Albedo, February-October (mean of the monthly means; a month needs >= %.0f %% of the ' ...
        'period''s years), and its difference from the %s normal; statewide values over the %d cells with every ' ...
        'period; grey: no value in that period'], ...
        100 * E.perMinFrac, bs, R.perRegN(10)), C);
    saveFig(f, fullfile(E.figDir, 'fig07_five_year_means.pdf'));
end

% ---------- fig09 (optional, E.snowYears): statewide daily Albedo and NDSI snow cover, 1 Sep - 30 Aug ----------
% one figure, the snow years stacked on the same date axis: every date is drawn at its own day and month of the
% reference snow year 1 Sep 2003 - 30 Aug 2004 (2004 has a 29 Feb, so no date moves); ticks on the 1st of each month
if any(E.figures == 9) && ~isempty(X.snow)
    nS = numel(X.snow);
    f = newFig(13, 0.8 + 5.6 * nS, E);
    t = tiles(f, 3 * nS, 1);
    ca = C.blue(6, :);  cs = C.orng(5, :);
    xl = datenum([2003 2004], [9 8], [1 30]);  tk = datenum(2003, 9:20, 1);
    tl = cellfun(@(mo) ['1 ' mo], cellstr(datestr(tk, 'mmm')), 'UniformOutput', false);
    for k = 1:nS
        sw = X.snow(k);
        dv = datevec(sw.dn);
        xs = datenum(2003 + dv(:, 1) - sw.y0, dv(:, 2), dv(:, 3));   % the same day and month in 2003-04
        dr = datevec(xs);
        assert(isequal(dr(:, 2:3), dv(:, 2:3)) && xs(1) == xl(1) && xs(end) == xl(2), ...
            'fig09: the dates of %d-%d do not map onto 1 Sep 2003 - 30 Aug 2004', sw.y0, sw.y0 + 1);
        a = tileAt(t, 3 * k - 2, [2 1], C);
        plot(a, xs, sw.albClim, '--', 'Color', 0.45 * ca + 0.55, 'LineWidth', 1.2);
        plot(a, xs, sw.ndsiClim, '--', 'Color', 0.45 * cs + 0.55, 'LineWidth', 1.2);
        plot(a, xs, sw.alb, '-', 'Color', 0.45 * ca + 0.55, 'LineWidth', 0.6);
        plot(a, xs, sw.ndsi, '-', 'Color', 0.45 * cs + 0.55, 'LineWidth', 0.6);
        h1 = plot(a, xs, sw.albS, '-', 'Color', ca, 'LineWidth', 2);
        h2 = plot(a, xs, sw.ndsiS, '-', 'Color', cs, 'LineWidth', 2);
        h3 = plot(a, NaN, NaN, '--', 'Color', C.ink2, 'LineWidth', 1.2);
        set(a, 'XLim', xl, 'YLim', [0 1], 'XTick', tk, 'XTickLabel', tl, 'XGrid', 'on', 'YGrid', 'on');
        ylabel(a, 'statewide mean (0-1)');
        if k == 1
            legend(a, [h1 h2 h3], {'blue-sky Albedo', 'NDSI snow cover', sprintf('mean of %s (same colour)', bs)}, ...
                'Box', 'off', 'Location', 'best', 'Orientation', 'horizontal', 'FontSize', 9);   % least overlap
        end
        panelTitle(a, sprintf('1 Sep %d - 30 Aug %d', sw.y0, sw.y0 + 1), C, sprintf(['thin: daily; thick: 7-day mean; ' ...
            'a day needs >= %.0f %% of the AOI pixels valid (Albedo: %d days, NDSI: %d days)'], 100 * E.minRegDayCov, ...
            nnz(isfinite(sw.alb)), nnz(isfinite(sw.ndsi))));
        % share of the AOI pixels behind each daily value: in the dark, low-sun months the 'statewide' mean rests
        % on the few (mostly southern, clear) pixels with a value, so it is not representative of the state there
        b = tileAt(t, 3 * k, [1 1], C);
        plot(b, xs, 100 * sw.albCover, '-', 'Color', ca, 'LineWidth', 1);
        plot(b, xs, 100 * sw.ndsiCover, '-', 'Color', cs, 'LineWidth', 1);
        plot(b, xl, 100 * E.minRegDayCov * [1 1], ':', 'Color', C.ink2, 'LineWidth', 1);
        set(b, 'XLim', xl, 'YLim', [0 100], 'XTick', tk, 'XTickLabel', tl, 'XGrid', 'on', 'YGrid', 'on');
        ylabel(b, 'AOI pixels valid (%)');
        panelTitle(b, sprintf('coverage %d-%d', sw.y0, sw.y0 + 1), C, sprintf(['share of the AOI pixels with a value ' ...
            '(dotted: the %.0f %% minimum)'], 100 * E.minRegDayCov));
    end
    sy = strjoin(arrayfun(@(y) sprintf('%d-%d', y, y + 1), [X.snow.y0], 'UniformOutput', false), ' and ');
    layoutTitle(t, sprintf(['Statewide daily blue-sky Albedo and MOD10A1 NDSI snow cover, 1 Sep - 30 Aug of %s ' ...
        '(pooled mean over the valid 1-km AOI pixels; the same date axis for every snow year)'], sy), C);
    saveFig(f, fullfile(E.figDir, 'fig09_snow_years.pdf'));
end
end

function drawCycleFig(cycYear, cycMean, met, ylim, ylab, what, src, fn, cyc, C, E)
% fig01 layout (statewide panel spanning 3 x 2 tiles, then the 9 ecoregions): every baseline year's 7-day
% smoothed regional curve and their mean; cycYear [10 x 365 x years], cycMean [10 x 365]. Without groups
% (cyc.grp empty) every year has its own colour (cyc.ycol) and the figure a colour bar. With groups
% (E.yearGroups, yearGroups) every year is a pale curve in its group's colour, every group gets its mean curve,
% and every panel's legend gives, per group and for all years, the mean over the years of the annual metrics
% met.show of met.val [10 x years x metrics] (groupLegend)
yb = cyc.yb;  nb = numel(yb);  names = cyc.names;  P = cyc.P;  xd = cyc.xd;  win = cyc.win;
monTick = cyc.monTick;  bs = cyc.bs;  grp = cyc.grp;
assert(size(cycYear, 3) == nb && isequal(size(cycMean), [10 365]), 'drawCycleFig: curves of the %d baseline years expected', nb);
if isempty(grp)
    ttl = sprintf('Annual cycle of daily %s, every year %s and their mean (%s)', what, bs, src);
else
    assert(isequal(size(met.val), [10 nb numel(met.name)]), 'drawCycleFig: metrics of the %d baseline years expected', nb);
    ttl = sprintf(['Annual cycle of daily %s %s (%s): every year (pale, in the colour of its group of years), the mean ' ...
        'of each group (colour; days with >= %.0f %% of its years) and the mean of all years (black)'], what, bs, src, ...
        100 * E.groupMinFrac);
    if ~isempty(met.show)
        ttl = sprintf(['%s. Legend: mean over the years of each year''s %s (from the year''s curve, computed as for ' ...
            'tab11; -: values in < %.0f %% of the years)'], ttl, metricText(met), 100 * E.groupMinFrac);
    end
end
f = newFig(30, 20, E);
t = tiles(f, 3, 5);
layoutTitle(t, ttl, C);
if ~isempty(grp), mono = monoFont(); end      % legend tables: columns aligned in a fixed-width font

for p = 1:10
    j = P(p);
    if p == 1
        a = tileAx(t, 0, [3 2], C);  span = [3 2];
    else
        a = tileAx(t, p, [1 1], C);  span = [1 1];
    end

    if isempty(grp)
        for i = 1:nb
            plot(a, xd(win), double(cycYear(j, win, i)), '-', 'Color', cyc.ycol(i, :), 'LineWidth', 0.7);
        end
    else
        for i = 1:nb                    % pale years: the edge of a patch (MATLAB lines have no documented
            y = double(cycYear(j, win, i));   % transparency); the NaN at the end keeps the curve open
            patch(a, [xd(win) NaN], [y NaN], 'k', 'FaceColor', 'none', 'EdgeColor', grp.col(grp.id(i), :), ...
                'EdgeAlpha', E.yearAlpha, 'LineWidth', 0.7);
        end
        nG = numel(grp.n);  hg = gobjects(1, nG);  GM = nan(nG, nnz(win));
        for g = 1:nG                    % group mean: mean of its years' curves, on days with >= E.groupMinFrac
            Yg = double(reshape(cycYear(j, :, grp.id == g), 365, grp.n(g)));
            gm = nanMean(Yg, 2)';  gm(sum(isfinite(Yg), 2)' < E.groupMinFrac * grp.n(g)) = NaN;
            hg(g) = plot(a, xd(win), gm(win), '-', 'Color', grp.col(g, :), 'LineWidth', 1.6);
            GM(g, :) = gm(win);
        end
    end
    plot(a, xd(win), cycMean(j, win), '-', 'Color', 'w', 'LineWidth', 4.2);   % white edge: the mean stays
    hl = plot(a, xd(win), cycMean(j, win), '-', 'Color', C.ink, 'LineWidth', 2.2);   % readable over the years
    set(a, 'XLim', [32 304], 'YLim', ylim, monTick{:}, 'XGrid', 'on', 'YGrid', 'on', 'FontSize',13);
    panelTitle(a, names{j}, C);
    if ~any([1 2 5 8]==p), set(a, 'YTickLabel', {}); end
    if p <= 7 && p~=1, set(a, 'XTickLabel', {}); end
    if isempty(grp)
        if p == 1, legend(a, hl, {sprintf('mean %s', bs)}, 'Box', 'off', 'Location', 'northeast', 'FontSize', E.cycleLegendFont); end
    else                                % every panel: a table of the groups' metrics (groupLegend)
        [hdr, lab] = groupLegend(met, j, grp, bs, p == 1, E.groupMinFrac);
        h = [hg hl];
        if ~isempty(hdr)                % header row: an entry without a symbol
            h = [plot(a, NaN, NaN, 'LineStyle', 'none', 'Color', 'w'), h];  lab = [{hdr}, lab];
        end
        Yp = [double(reshape(cycYear(j, win, :), nnz(win), nb))'; GM; cycMean(j, win)];   % every curve drawn
        wh = legendShare(numel(lab), max(cellfun(@numel, lab)), E.cycleLegendFont, f, [3 5], span);
        if any(wh > 1)
            warning('alb_eda:legendSize', ['fig01 %s: the legend may not fit its panel; lower E.cycleLegendFont ' ...
                'or enlarge the figure (newFig in drawCycleFig)'], names{j});
        end
        lg = legend(a, h, strrep(lab, ' ', char(160)), 'Location', emptyCorner(xd(win), Yp, [32 304], ylim, wh), ...
            'FontName', mono, 'FontSize', E.cycleLegendFont, 'Interpreter', 'latex', 'Box', 'on', 'Color', 'w', ...
            'EdgeColor', C.grid);       % no-break spaces: the padding of the columns is kept as typed
        lg.ItemTokenSize = [18 18];     % shorter symbols (default 30 points)
    end
end
ylabel(t,ylab,'Interpreter','latex', 'FontSize',15);
xlabel(t,'','Interpreter','latex');
linkaxes(findall(t,'Type','axes'),'xy');

if isempty(grp)                         % one band per year, labelled at its centre
    layoutColorbar(a, cyc.ycol, [yb(1) - 0.5, yb(end) + 0.5], '', 'east', C, yb);
end

saveFig(f, fn);
end

function drawCycleBoth(cycA, meanA, cycN, meanN, hmA, hmN, fn, cyc, C, E)
% fig16 (E.figures 16): the fig01 layout (statewide panel spanning 3 x 2 tiles, then the 9 ecoregions) with blue-sky
% Albedo (solid) and NDSI snow cover (dashed) on one axis, both 0-1 (E.bothYLim, no yyaxis): the mean curve of each
% group of years (cyc.grp; days with >= E.groupMinFrac of its years) and of all years (black, white edge); the pale
% single years only with E.bothYears. cycA, cycN [10 x 365 x years], meanA, meanN [10 x 365]. Every legend (as in fig01,
% without a legend title): the mean over the years of the half-melt date of each year's Albedo (hmA) and NDSI (hmN)
% curve [10 x years], per group and for all years ('-' with values in < E.groupMinFrac of the years)
yb = cyc.yb;  nb = numel(yb);  names = cyc.names;  P = cyc.P;  xd = cyc.xd;  win = cyc.win;  grp = cyc.grp;
nG = numel(grp.n);  yl = E.bothYLim;
assert(isequal(size(cycA), [10 365 nb]) && isequal(size(cycN), [10 365 nb]) && isequal(size(meanA), [10 365]) && ...
    isequal(size(meanN), [10 365]) && isequal(size(hmA), [10 nb]) && isequal(size(hmN), [10 nb]), ...
    'drawCycleBoth: curves and half-melt dates of the %d baseline years expected', nb);
ttl = sprintf(['Annual cycle of daily blue-sky Albedo (solid) and MOD10A1 NDSI snow cover (dashed), both 0-1, %s ' ...
    '(regional daily means, 7-day smoothed): the mean of each group of years (colour; days with >= %.0f %% of its ' ...
    'years) and of all years (black). Legend: mean over the years of the half-melt date of each year''s Albedo and ' ...
    'NDSI curve (-: values in < %.0f %% of the years)'], cyc.bs, 100 * E.groupMinFrac, 100 * E.groupMinFrac);
f = newFig(30, 20, E);
t = tiles(f, 3, 5);
layoutTitle(t, ttl, C);
mono = monoFont();
right = @(c, n) [repmat(' ', 1, n - numel(c)), c];       % c right-aligned in n characters
ls = {'-', '--'};  Ys = {cycA, cycN};  Ms = {meanA, meanN};
for p = 1:10
    j = P(p);
    if p == 1
        a = tileAx(t, 0, [3 2], C);  span = [3 2];
    else
        a = tileAx(t, p, [1 1], C);  span = [1 1];
    end
    hg = gobjects(1, nG);  Yp = [];
    for s = 1:2                                           % s = 1 Albedo (solid), 2 NDSI snow cover (dashed)
        Y = Ys{s};
        if E.bothYears
            for i = 1:nb
                y = double(Y(j, win, i));
                patch(a, [xd(win) NaN], [y NaN], 'k', 'FaceColor', 'none', 'EdgeColor', grp.col(grp.id(i), :), ...
                    'EdgeAlpha', E.yearAlpha, 'LineWidth', 0.6, 'LineStyle', ls{s});
            end
        end
        for g = 1:nG
            Yg = double(reshape(Y(j, :, grp.id == g), 365, grp.n(g)));
            gm = nanMean(Yg, 2)';  gm(sum(isfinite(Yg), 2)' < E.groupMinFrac * grp.n(g)) = NaN;
            h = plot(a, xd(win), gm(win), ls{s}, 'Color', grp.col(g, :), 'LineWidth', 1.6, 'Marker', 'none');
            if s == 1, hg(g) = h; end
            Yp = [Yp; gm(win)];                                                       %#ok<AGROW>
        end
        mm = Ms{s}(j, win);
        plot(a, xd(win), mm, '-', 'Color', 'w', 'LineWidth', 4.2, 'Marker', 'none');   % white edge
        h = plot(a, xd(win), mm, ls{s}, 'Color', C.ink, 'LineWidth', 2.2, 'Marker', 'none');
        if s == 1, hl = h; end
        Yp = [Yp; mm];                                                                %#ok<AGROW>
    end
    set(a, 'XLim', [32 304], 'YLim', yl, cyc.monTick{:}, 'XGrid', 'on', 'YGrid', 'on');
    panelTitle(a, names{j}, C);
    if ~any([1 2 5 8] == p), set(a, 'YTickLabel', {}); end           % ticks: first column
    if p <= 7 && p ~= 1, set(a, 'XTickLabel', {}); end
    if p == 1, ylabel(a, 'Albedo (solid), NDSI snow cover (dashed), 0-1'); end
    % legend table: per group and for all years, the mean half-melt date of Albedo and of NDSI
    nm = [grp.label, {cyc.bs}];  val = cell(nG + 1, 2);
    for g = 1:nG + 1
        if g <= nG, in = grp.id == g; else, in = true(1, nb); end
        if p == 1, nm{g} = sprintf('%s (%d y)', nm{g}, nnz(in)); end
        H = {hmA(j, in), hmN(j, in)};
        for k = 1:2
            vv = H{k}(isfinite(H{k}));
            if numel(vv) >= E.groupMinFrac * nnz(in), val{g, k} = doyLabel(mean(vv)); else, val{g, k} = '-'; end
        end
    end
    w = max(cellfun(@numel, [val(:)', {'Albedo'}]));
    hdr = [right('Albedo', w), '  ', right('NDSI', w)];
    lab = cell(1, nG + 1);
    for g = 1:nG + 1, lab{g} = [right(val{g, 1}, w), '  ', right(val{g, 2}, w)]; end
    if p == 1                                               % the statewide panel names the years of each row
        wl = max(cellfun(@numel, nm));
        hdr = [repmat(' ', 1, wl + 2), hdr];
        for g = 1:nG + 1, lab{g} = [nm{g}, repmat(' ', 1, wl - numel(nm{g}) + 2), lab{g}]; end
    end
    h = [plot(a, NaN, NaN, 'LineStyle', 'none', 'Color', 'w'), hg, hl];  lab = [{hdr}, lab];
    wh = legendShare(numel(lab), max(cellfun(@numel, lab)), E.cycleLegendFont, f, [3 5], span);
    if any(wh > 1)
        warning('alb_eda:legendSize', ['fig16 %s: the legend may not fit its panel; lower E.cycleLegendFont ' ...
            'or enlarge the figure (newFig in drawCycleBoth)'], names{j});
    end
    lg = legend(a, h, strrep(lab, ' ', char(160)), 'Location', emptyCorner(xd(win), Yp, [32 304], yl, wh), ...
        'FontName', mono, 'FontSize', E.cycleLegendFont, 'Interpreter', 'latex', 'Box', 'on', 'Color', 'w', ...
        'EdgeColor', C.grid);       % as fig01: no-break spaces keep the padding of the columns
    lg.ItemTokenSize = [18 18];
end
linkaxes(findall(t, 'Type', 'axes'), 'x');
saveFig(f, fn);
end

function drawSeriesFig(cycA, cycD, yb, names, P, fn, C, E)
% fig17 (E.figures 17): the daily regional means (7-day smoothed) of every baseline year joined into one series per
% region: blue-sky Albedo, NDSI snow cover and SWE scaled by its largest value in the region (left axis, 0-1) and
% 2-m air temperature (right axis, E.seriesTLim). cycA [10 x 365 x years]; cycD [10 x 365 x years x 3] = NDSI, t2m
% (degC), SWE (mm). Panel order P: the statewide panel first (two rows high), then the nine ecoregions
nb = numel(yb);
assert(size(cycA, 2) == 365 && size(cycA, 3) == nb && isequal(size(cycD, 2:4), [365 nb 3]), ...
    'drawSeriesFig: curves [10 x 365 x %d] expected', nb);
tt = reshape(yb(:)' + ((1:365)' - 0.5) / 365, 1, []);   % decimal year of every day (365-day years), year after year
cAlb = C.ink;  cNdsi = [0.106 0.604 0.667];  cT = [0.851 0.345 0.169];  cSwe = [0.561 0.663 0.761];
f = newFig(24, 26, E);
nP = numel(P);
t = tiles(f, nP + 1, 1);
layoutTitle(t, sprintf(['Daily regional means %d-%d (7-day smoothed): blue-sky Albedo, MOD10A1 NDSI snow cover and ' ...
    'ERA5-Land snow water equivalent scaled by its largest value in the region (left axis), and ERA5-Land 2-m air ' ...
    'temperature (right axis). No Albedo or NDSI in the polar night (November-January)'], yb(1), yb(end)), C);
dg = char(176);
for p = 1:nP
    j = P(p);
    if p == 1, a = tileAt(t, 1, [2 1], C); else, a = tileAt(t, p + 1, [1 1], C); end
    alb = reshape(double(cycA(j, :, :)), 1, []);
    nd = reshape(double(cycD(j, :, :, 1)), 1, []);
    tk = reshape(double(cycD(j, :, :, 2)), 1, []);
    sw = reshape(double(cycD(j, :, :, 3)), 1, []);
    smax = max(sw(isfinite(sw)));
    yyaxis(a, 'left');  hold(a, 'on');
    if ~isempty(smax) && smax > 0
        s = sw / smax;  s(~isfinite(s)) = 0;
        hS = area(a, tt, s, 'FaceColor', cSwe, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
        sub = sprintf('SWE scaled by its largest value, %.0f mm', smax);
    else
        hS = plot(a, NaN, NaN, '-', 'Color', cSwe);  sub = 'no SWE';
    end
    hN = plot(a, tt, nd, '-', 'Color', cNdsi, 'LineWidth', 0.9, 'Marker', 'none');
    hA = plot(a, tt, alb, '-', 'Color', cAlb, 'LineWidth', 1.0, 'Marker', 'none');
    set(a, 'YLim', [0 1]);
    a.YAxis(1).Color = C.ink2;
    yyaxis(a, 'right');  hold(a, 'on');
    plot(a, [yb(1) yb(end) + 1], [0 0], ':', 'Color', cT, 'LineWidth', 0.7, 'Marker', 'none');
    hT = plot(a, tt, tk, '-', 'Color', cT, 'LineWidth', 0.6, 'Marker', 'none');
    set(a, 'YLim', E.seriesTLim);
    a.YAxis(2).Color = cT;
    yyaxis(a, 'left');
    xt = yb(1):yb(end) + 1;
    set(a, 'XLim', [yb(1) yb(end) + 1], 'XTick', xt, 'XGrid', 'on', 'YGrid', 'on');
    if p == 1 || p == nP
        set(a, 'XTickLabel', [arrayfun(@num2str, yb(:)', 'UniformOutput', false), {''}]);
    else
        set(a, 'XTickLabel', {});
    end
    panelTitle(a, names{j}, C, sub);
    if p == 1
        ylabel(a, 'Albedo, NDSI, SWE / max');
        yyaxis(a, 'right');  ylabel(a, ['air temperature (' dg 'C)']);  yyaxis(a, 'left');
        legend(a, [hA hN hS hT], {'blue-sky Albedo (left axis)', 'NDSI snow cover, 0-1 (left axis)', ...
            'SWE / its largest value in the region (left axis)', ['2-m air temperature, ' dg 'C (right axis)']}, ...
            'Orientation', 'horizontal', 'Location', 'northoutside', 'FontSize', 9, 'Box', 'off');
    end
end
linkaxes(findall(t, 'Type', 'axes'), 'x');
saveFig(f, fn);
end

function yl = niceLimits(v)
% axis limits around the values v on a 1-2-5 step (about 5 steps)
lo = min(v);  hi = max(v);
if hi <= lo, yl = lo + [-1 1]; return; end
st = 10 ^ floor(log10((hi - lo) / 5));
for m = [1 2 5 10], if (hi - lo) / (m * st) <= 6, st = m * st; break; end, end
yl = [floor(lo / st) * st, ceil(hi / st) * st];
end

function grp = yearGroups(Gy, yb, anchors)
% groups of the baseline years yb in the fig01 annual cycles. Gy (E.yearGroups): one row [first last] per group,
% in order and together covering every baseline year. -> grp.id (group of every year), grp.n (years per group),
% grp.label (the first and last baseline year of the group) and grp.col (the anchor colours of E.yearPalette
% spread evenly over the groups: as many groups as anchors = the anchors themselves)
nG = size(Gy, 1);
assert(nG >= 1 && size(Gy, 2) == 2 && all(Gy(:, 1) <= Gy(:, 2)) && all(Gy(2:end, 1) > Gy(1:end - 1, 2)), ...
    'E.yearGroups: one row [first last] per group, in order and not overlapping');
grp.id = zeros(1, numel(yb));
for g = 1:nG, grp.id(yb >= Gy(g, 1) & yb <= Gy(g, 2)) = g; end
assert(all(grp.id > 0), 'E.yearGroups leaves out the baseline years %s', mat2str(yb(grp.id == 0)));
grp.n = accumarray(grp.id(:), 1, [nG 1])';
assert(all(grp.n > 0), 'E.yearGroups: a group without baseline years');
grp.label = cell(1, nG);
for g = 1:nG
    y = yb(grp.id == g);
    if numel(y) == 1, grp.label{g} = sprintf('%d', y); else, grp.label{g} = sprintf('%d-%d', y(1), y(end)); end
end
grp.col = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, nG));
end

function q = groupMetrics(E, v)
% numbers (cycleMetrics) of the annual metrics in the fig01 legends of variable v (E.groupMetrics)
k = find(strcmp(E.groupMetrics(:, 1), v));
assert(numel(k) == 1, 'E.groupMetrics needs one row for ''%s''', v);
q = E.groupMetrics{k, 2};  q = q(:)';
assert(all(ismember(q, 1:10)), 'E.groupMetrics for ''%s'': numbers 1-10 of cycleMetrics', v);
end

function [hdr, lab] = groupLegend(met, j, grp, bs, long, minFrac)
% fig01 legend table of region j, in columns for a fixed-width font: the header hdr (short names of the metrics
% met.show) and one row per group and one for all years (lab) with the mean over their years of each metric
% ('-' with values in < minFrac of the years), right-aligned. long (statewide panel, the key of the figure):
% every row starts with the years of its group and their number, '2003-2009 (7 y)'. Without metrics, hdr is
% '' and the rows are the years (long: with their number).
nG = numel(grp.n);  nk = numel(met.show);  nm = cell(1, nG + 1);  val = cell(nG + 1, nk);
for g = 1:nG + 1
    if g <= nG, in = grp.id == g;  nm{g} = grp.label{g}; else, in = true(size(grp.id));  nm{g} = bs; end
    if long, nm{g} = sprintf('%s (%d y)', nm{g}, nnz(in)); end
    for k = 1:nk
        v = met.val(j, in, met.show(k));  v = v(isfinite(v));
        if numel(v) >= minFrac * nnz(in), val{g, k} = metricValue(mean(v), met.unit{met.show(k)}, met.fmt); else, val{g, k} = '-'; end
    end
end
hdr = '';  lab = nm;
if nk == 0, return; end
head = metricShort(met.show, met.name);
w = zeros(1, nk);
for k = 1:nk, w(k) = max(cellfun(@numel, [val(:, k)', head(k)])); end
right = @(c, n) [repmat(' ', 1, n - numel(c)), c];      % c right-aligned in n characters
hdr = strjoin(arrayfun(@(k) right(head{k}, w(k)), 1:nk, 'UniformOutput', false), '  ');
for g = 1:nG + 1
    lab{g} = strjoin(arrayfun(@(k) right(val{g, k}, w(k)), 1:nk, 'UniformOutput', false), '  ');
end
if long                                                    % the years of every row in front
    wl = max(cellfun(@numel, nm));
    hdr = [repmat(' ', 1, wl + 2), hdr];
    for g = 1:nG + 1, lab{g} = [nm{g}, repmat(' ', 1, wl - numel(nm{g}) + 2), lab{g}]; end
end
end

function s = metricShort(q, name)
% column names of the cycleMetrics metrics q in the fig01 legend tables (name: cycleMetrics' names, checked)
short = {'Mar-Oct', 'peak', 'peak date', 'JJA', 'JJA min', 'amplitude', 'half-melt', 'melt days', 'return', 'snow-free'};
assert(numel(name) == numel(short) && strcmp(name{7}, 'half-melt date'), 'metricShort: cycleMetrics has changed');
s = short(q);
end

function fn = monoFont()
% fixed-width font of the fig01 legend tables: Consolas where installed (Windows), else MATLAB's own
if any(strcmpi(listfonts, 'Consolas')), fn = 'Consolas'; else, fn = get(groot, 'FixedWidthFontName'); end
end

function wh = legendShare(nr, nc, fs, f, grid, span)
% approximate size of a legend of nr rows of nc fixed-width characters (font fs pt) as a share [width height] of
% an axes spanning span = [rows cols] tiles of a grid = [rows cols] layout filling the figure f: symbol, gaps
% and characters (0.6 em wide, rows 1.3 em) over the axes size (taken as 95 % of the figure width and 88 % of
% its height, shared by the tiles); > 1: larger than the axes
u = get(f, 'Units');  set(f, 'Units', 'inches');  fp = get(f, 'Position');  set(f, 'Units', u);
w = 0.95 * fp(3) * span(2) / grid(2);  h = 0.88 * fp(4) * span(1) / grid(1);      % inches
wh = [(18 + 24 + 0.6 * fs * nc) / 72 / w, (1.3 * fs * nr + 8) / 72 / h];
end

function loc = emptyCorner(x, Y, xl, yl, wh)
% legend location: the corner of the axes (limits xl, yl) whose box of wh = [width height] (share of the axes)
% plus a margin of 3 % holds the fewest visible points of the curves Y (rows, at x); ties go to the first of
% NE, NW, SE, SW
wh = min(1, wh + 0.03);
xn = repmat((x(:)' - xl(1)) / diff(xl), size(Y, 1), 1);  yn = (Y - yl(1)) / diff(yl);
ok = isfinite(yn) & yn >= 0 & yn <= 1 & xn >= 0 & xn <= 1;  xn = xn(ok);  yn = yn(ok);
right = xn >= 1 - wh(1);  left = xn <= wh(1);  top = yn >= 1 - wh(2);  bottom = yn <= wh(2);
cnt = [nnz(right & top), nnz(left & top), nnz(right & bottom), nnz(left & bottom)];
locs = {'northeast', 'northwest', 'southeast', 'southwest'};
[~, k] = min(cnt);  loc = locs{k};
end

function s = metricValue(v, unit, fmt)
% a metric in a fig01 legend: a date ('13 May'), a number of days ('96 d') or a value (format fmt)
switch unit
    case 'day of year', s = doyLabel(v);
    case 'days',        s = sprintf('%.0f d', v);
    otherwise,          s = sprintf(fmt, v);
end
end

function s = metricText(met)
% the legend metrics met.show for the fig01 titles: 'mean Mar-Oct and half-melt date' (+ the unit of values)
p = met.name(met.show);
for k = 1:numel(p)
    if ~isempty(met.vunit) && ~any(strcmp(met.unit{met.show(k)}, {'day of year', 'days'}))
        p{k} = sprintf('%s (%s)', p{k}, met.vunit);
    end
end
if numel(p) <= 2, s = strjoin(p, ' and '); else, s = [strjoin(p(1:end - 1), ', ') ' and ' p{end}]; end
end

function s = doyLabel(d)
% day of year (365-day calendar) -> '2 May'
s = strtrim(datestr(datenum(2001, 1, 0) + round(d), 'dd mmm'));
if s(1) == '0', s = s(2:end); end
end

function A = yearPalette(name)
% anchor colours (rows: first to last) of the fig01 annual cycles; drawFigures spreads the groups of
% E.yearGroups (yearGroups), or the years (cmapFrom), evenly over them. The sequential ones run from OKLab
% lightness 0.76 (2:1 contrast on white: a thin line stays visible) to 0.37 (clear of the black mean curve),
% in even steps.
switch name
    case 'tealIndigo'                   % OKLCH hue 180 -> 275, chroma 0.10
        h = {'#5fc6b3', '#4bb9b2', '#39abb0', '#2b9ead', '#278fa8', '#2281a3', '#23729d', '#286495', '#2c558b', ...
             '#304780', '#313974'};
    case 'viridis'                      % matplotlib viridis from 0.761 to 0.121 of the map
        h = {'#63cb5f', '#42be71', '#2ab07f', '#1fa188', '#20928c', '#26828e', '#2c728e', '#33638d', '#3b518b', ...
             '#423f85', '#472a7a'};
    case 'batlow'                       % Crameri batlow (Crameri et al. 2020) from 0.769 to 0.132 of the map
        h = {'#f6a077', '#e39754', '#c8913b', '#a78b2c', '#88842f', '#6c7c3c', '#52744a', '#386c58', '#246160', ...
             '#185462', '#114460'};
    case 'blueRed'                      % green, blue, red (one per decade of E.yearGroups). The green was
        % [0.73 1.00 0.40] (#BAFF66): 1.2:1 contrast on white, its pale year curves vanish. #22984C (OKLCH
        % lightness 0.60) has 3.7:1 and stays apart from the red for red-green colour-blind readers
        A = [0.133, 0.596, 0.298; 0, 0.3, 0.99; 1.0, 0.19, 0.19];  return
    otherwise
        error('E.yearPalette = ''%s'': use ''tealIndigo'', ''viridis'', ''batlow'' or ''blueRed''', name);
end
A = reshape(sscanf(strjoin(cellfun(@(c) c(2:end), h, 'UniformOutput', false), ''), '%2x'), 3, [])' / 255;
end

function C = palette()
hex = @(s) reshape(sscanf(strjoin(cellfun(@(c) c(2:end), s, 'UniformOutput', false), ''), '%2x'), 3, [])' / 255;
C.ink = hex({'#1f1f1d'});  C.ink2 = hex({'#52514e'});  C.grid = hex({'#e4e2dc'});
C.never = hex({'#b5b3ac'});
C.blue = hex({'#eef5fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'});
C.orng = hex({'#fef3ec', '#fde5d8', '#f9c0a3', '#f39a70', '#eb6834', '#c94f1f', '#9c3a14', '#6e280d'});
C.div  = hex({'#7d1a1a', '#b52d2c', '#e34948', '#f3a19b', '#fbd5d2', '#f0efec', '#cde2fb', '#86b6ef', '#3987e5', ...
              '#1c5cab', '#0d366b'});                  % red = darker, grey midpoint, blue = brighter
C.seq = cmapFrom(C.blue);  C.divmap = cmapFrom(C.div);
end

function M = cmapFrom(anchors)
t = linspace(0, 1, size(anchors, 1));
M = interp1(t, anchors, linspace(0, 1, 256));
end

function c = textOn(dark, C)
% text colour on a heatmap cell
if dark, c = [1 1 1]; else, c = C.ink; end
end

function f = newFig(w, h, E)
f = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 w h], 'Visible', E.figVisible, ...
    'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
if isprop(f, 'Theme'), f.Theme = 'light'; end         % R2025a+: never inherit a dark theme
end

function styleAx(ax, C)
set(ax, 'Box', 'off', 'TickDir', 'in', 'FontSize', 9, 'XColor', C.ink2, 'YColor', C.ink2, 'Layer', 'top', ...
    'GridColor', C.grid, 'GridAlpha', 1, 'Color', 'w');
hold(ax, 'on');
end

function t = tiles(f, nr, nc)
% tiled chart layout (nr x nc tiles) filling the figure
t = tiledlayout(f, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function t = tilesIn(t0, nr, nc, tile, span)
% nr x nc layout nested in tile 'tile' of the layout t0, spanning span = [rows cols] of its tiles
t = tiledlayout(t0, nr, nc, 'TileSpacing', 'compact', 'Padding', 'none');
t.Layout.Tile = tile;  t.Layout.TileSpan = span;
end

function ax = tileAx(t, tile, span, C)
% styled axes in tile 'tile' of the layout t (tiles numbered row by row), spanning span = [rows cols]
if tile == 0
    ax = nexttile(t, span);
else
    ax = nexttile(t);
end
styleAx(ax, C);
end

function ax = tileAt(t, tile, span, C)
% styled axes at tile number 'tile' (row by row) of the layout t, spanning span = [rows cols]
ax = nexttile(t, tile, span);
styleAx(ax, C);
end

function layoutTitle(t, str, C)
%LAYOUTTITLE - Set a wrapped bold title for an entire layout
%
% Input arguments:
%   t   - tiled layout handle
%   str - title text to wrap
%   C   - palette struct with ink color
%
% Output arguments:
%   none
% title of a whole layout, wrapped to the figure width (a layout title does not wrap by itself)
% Measure figure width in inches so wrapText can estimate line length.
f = ancestor(t, 'figure');  u = f.Units;  f.Units = 'inches';  w = f.Position(3);  f.Units = u;
% Wrap the title to the figure width before applying layout title formatting.
h = title(t, wrapText(str, floor(w * 72 / (0.6 * 12))));
% Use bold, 12-point ink-colored text for the final layout title.
set(h, 'FontWeight', 'bold', 'FontSize', 12, 'Color', C.ink);
end

function lines = wrapText(str, n)
% str split at its line breaks, each line wrapped at spaces to <= n characters -> cell array of lines
parts = regexp(str, '\n', 'split');  lines = {};
for p = 1:numel(parts)
    w = strsplit(parts{p}, ' ');  cur = '';
    for k = 1:numel(w)
        if isempty(cur), cand = w{k}; else, cand = [cur ' ' w{k}]; end
        if numel(cand) > n && ~isempty(cur)
            lines{end + 1} = cur;  cur = w{k};                                         %#ok<AGROW>
        else
            cur = cand;
        end
    end
    lines{end + 1} = cur;                                                              %#ok<AGROW>
end
end

function panelTitle(ax, str, C, sub)
% Format the title and default the subtitle when omitted.
if nargin < 4, sub = ''; end
% Create the main and subtitle text objects with LaTeX-safe strings.
[ht, hs] = title(ax, latexSafe(str), latexSafe(sub));          % LaTeX: % and >= print as typed
set(ht, 'FontWeight', 'bold', 'FontSize', 12, 'Color', C.ink, 'Visible', 'on', 'Interpreter','latex');
set(hs, 'FontWeight', 'normal', 'FontSize', 10, 'Color', C.ink2, 'Visible', 'on', 'Interpreter','latex');
ax.TitleHorizontalAlignment = 'center';
set(ax, 'Box','on')
end

function s = latexSafe(s)
% text for the LaTeX interpreter: outside $...$, % & # _ are escaped (unless they already are), >= <=
% < > become math symbols, and the Greek letters \alpha ... \omega and the math-only commands (arrows, \pm, \times,
% \cdot, \circ, \approx, \sim, \infty, \leq, \geq, ...) are set in math mode; the parts inside $...$ are kept. With
% the LaTeX interpreter a bare % starts a comment (the rest of the text is lost), < > print as inverted marks and
% \beta or \rightarrow outside math mode is an error (no text). A cell array: every line
if iscell(s), s = cellfun(@latexSafe, s, 'UniformOutput', false); return; end
if isempty(s), return; end
assert(ischar(s), 'latexSafe: text expected');
p = regexp(s, '\$', 'split');
assert(mod(numel(p), 2) == 1, 'latexSafe: unbalanced $ in ''%s''', s);
gk = ['\\(alpha|beta|gamma|delta|epsilon|kappa|lambda|mu|nu|rho|sigma|tau|phi|xi|chi|psi|omega|Gamma|Delta|' ...
    'Sigma|Phi|Psi|Omega|rightarrow|leftarrow|leftrightarrow|Rightarrow|Leftarrow|to|uparrow|downarrow|pm|mp|' ...
    'times|cdot|div|circ|approx|sim|simeq|propto|infty|leq|geq|neq|ll|gg|partial|nabla)(?![A-Za-z])'];
for k = 1:2:numel(p)                                       % the parts outside math mode
    t = p{k};
    [u, m] = regexp(t, '(?<!\\)[%&#_]', 'split', 'match');
    t = u{1};  for i = 1:numel(m), t = [t, '\', m{i}, u{i + 1}]; end                          %#ok<AGROW>
    [u, m] = regexp(t, gk, 'split', 'match');
    t = u{1};  for i = 1:numel(m), t = [t, '$', m{i}, '$', u{i + 1}]; end                     %#ok<AGROW>
    t = strrep(strrep(t, '>=', '$\geq$'), '<=', '$\leq$');
    [u, m] = regexp(t, '[<>]', 'split', 'match');
    t = u{1};  for i = 1:numel(m), t = [t, '$', m{i}, '$', u{i + 1}]; end                     %#ok<AGROW>
    p{k} = t;
end
s = strjoin(p, '$');
end

function rgb = mapNoData(v, G, cmap, cl, C)
%MAPNODATA - Map kept-cell values and color missing data
%
% Input arguments:
%   v    - values for kept cells
%   G    - grid structure with keep mask and dimensions
%   cmap - colormap used for value mapping
%   cl   - color limits for mapping
%   C    - palette struct with reserved missing-data color
%
% Output arguments:
%   rgb - RGB image with missing cells colored separately
% map image of the values v of the kept cells through cmap over the limits cl (toRGB): kept cells without a value
% grey (C.never), the cells outside the AOI white
% Map kept-cell values, then mark missing values with the reserved grey.
rgb = toRGB(toGridC(v, G), cmap, cl, [1 1 1]);
% Build a mask for kept cells whose values are not finite.
gap = false(G.HC, G.WC);  gap(G.keep) = ~isfinite(v(:));
% Replace each RGB channel at missing cells with the "never observed" color.
for ch = 1:3, L = rgb(:, :, ch);  L(gap) = C.never(ch);  rgb(:, :, ch) = L; end
end

function cb = layoutColorbar(ax, cmap, cl, lab, where, C, ticks)
%LAYOUTCOLORBAR - Add a styled colorbar to a tiled layout or axes
%
% Input arguments:
%   ax    - target axes
%   cmap  - colormap to apply
%   cl    - color limits for the axes
%   lab   - colorbar label text
%   where - tile placement or 'axes' to keep it beside ax
%   C     - palette struct with ink colors
%   ticks - optional tick locations
%
% Output arguments:
%   cb - colorbar handle
% Keep axes colormap and limits aligned with the colorbar
colormap(ax, cmap);  set(ax, 'CLim', cl);
% Use caller-provided ticks when available
if nargin < 7, cb = colorbar(ax); else, cb = colorbar(ax, 'Ticks', ticks); end
% Reposition into the requested layout tile unless it stays beside the axes
if ~strcmp(where, 'axes'), cb.Layout.Tile = where; end
cb.Label.String = lab;
cb.Label.Interpreter = "latex";
set(cb, 'FontSize', 12, 'Color', C.ink2, 'TickDirection', 'out', 'TickLabelInterpreter','latex');
end

function layoutLegend(lg, where)
% move a legend to one side of its layout ('north', 'east', ...)
lg.Layout.Tile = where;
end

function rgb = toRGB(A, cmap, cl, nanColor)
%TORGB - Map values to RGB using a colormap and missing-value color
%
% Input arguments:
%   A        - numeric array of values to convert
%   cmap     - colormap with n rows and 3 columns
%   cl       - two-element color limits [low high]
%   nanColor - RGB triplet used for non-finite values
%
% Output arguments:
%   rgb - A reshaped RGB image array
% values -> RGB image through cmap and colour limits cl; NaN -> nanColor
% Determine colormap size and convert values to color indices.
n = size(cmap, 1);
k = round((A - cl(1)) / (cl(2) - cl(1)) * (n - 1)) + 1;
k = min(max(k, 1), n);
% Flag finite values and prefill the output with the missing-value color.
ok = isfinite(A);
rgb = zeros([numel(A) 3], 'single');
% Fill each RGB channel only where the input values are finite.
for ch = 1:3
    v = repmat(single(nanColor(ch)), numel(A), 1);
    v(ok(:)) = cmap(k(ok), ch);
    rgb(:, ch) = v;
end
% Restore the original array shape with a trailing RGB dimension.
rgb = reshape(rgb, [size(A) 3]);
end

function mapRGB(ax, rgb, x, y, G, C, cities)
%MAPRGB - Draw RGB map layers, boundaries, and optional city markers
%
% Input arguments:
%   ax     - target axes
%   rgb    - RGB image to display
%   x, y   - image coordinates
%   G      - grid metadata and contour fields
%   C      - color palette struct
%   cities - logical flag to draw city labels
image(ax, x, y, rgb);
set(ax, 'YDir', 'normal');  hold(ax, 'on');
% Overlay ecoregion contours
for r = 1:9
    contour(ax, G.xc, G.yc, double(G.ecoC == r), [0.5 0.5], 'LineColor', [0.44 0.43 0.41], 'LineWidth', 0.35);
end
% Draw the land-area boundary
contour(ax, G.xc, G.yc, double(G.lfC >= G.landFrac - 1e-9), [0.5 0.5], 'LineColor', [0.33 0.33 0.31], 'LineWidth', 0.5);
if cities
    % Mark each city with a point and label
    for c = 1:size(G.cities, 1)
        plot(ax, G.cities{c, 2} / 1e3, G.cities{c, 3} / 1e3, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', C.ink, 'MarkerSize', 4);
        text(ax, G.cities{c, 2} / 1e3 + 25, G.cities{c, 3} / 1e3 + 15, G.cities{c, 1}, 'FontSize', 7.5, 'Color', C.ink, ...
            'VerticalAlignment', 'bottom', 'Interpreter','latex');
    end
end
% Fix aspect ratio and hide axes decorations
set(ax, 'DataAspectRatio', [1 1 1], 'XLim', [G.x1(1) G.x1(end)], 'YLim', [G.y1(end) G.y1(1)], 'Visible', 'off');
end

function A = toGridC(v, G)
% values of the kept cells -> HC x WC map (NaN elsewhere)
A = nan(G.HC, G.WC);  A(G.keep) = v;
end

function saveFig(f, fn)
if exist('exportgraphics') > 0                                                  %#ok<EXIST>
    exportgraphics(f, fn, 'Resolution', 300);
else
    print(f, fn, '-dpng', '-r150');
end
close(f);
fprintf('  wrote %s\n', fn);
end

% ----------------------------- tables -----------------------------

function writeTables(R, G, E)
mn = R.mname;  nM = numel(mn);  regs = [G.order 0];
short = [G.short(G.order), {'All Alaska'}];
full = [G.regNames(G.order), {'All Alaska'}];
rg = [mn{1} mn{end}];  T = E.tabDir;
% tab01: dataset by year
yt = R.yearTable;  rows = cell(numel(yt), 9);
for i = 1:numel(yt)
    rows(i, :) = {yt(i).year, yt(i).bands, yt(i).first, yt(i).last, yt(i).valid_pixel_days, yt(i).valid_pct, ...
                  yt(i).magnitude_share_pct, yt(i).qa_mismatch, yt(i).in_baseline};
end
writeCsv(fullfile(T, 'tab01_dataset_by_year.csv'), {'year', 'bands', 'first', 'last', 'valid_pixel_days', ...
    ['valid_' rg '_pct'], ['magnitude_share_' rg '_pct'], 'qa_mismatch', 'in_baseline'}, rows);
% tab02: regional monthly climatology
rows = cell(10, 3 + nM);
for j = 1:10
    rows(j, :) = [{full{j}, R.regArea(j), R.regNever(j)}, num2cell(R.regClim(j, :))];
end
writeCsv(fullfile(T, 'tab02_regional_monthly_climatology.csv'), [{'region', 'area_km2', 'never_observed_pct'}, mn], rows);
% tab03: driver anomalies by region; tab03b: per cell (row-major order of the cell grid, or GeoTIFF maps)
rows = {};
for si = 1:numel(R.blkSeasons)
    for di = 1:numel(R.blkDrivers)
        r = squeeze(R.blkRho(si, di, :));  x = squeeze(R.blkXi(si, di, :));
        for j = 1:10
            if regs(j) > 0, sel = G.regC == regs(j); else, sel = G.regC > 0; end
            sel = sel & isfinite(r);
            if any(sel)
                v = {nnz(sel), median(r(sel)), pctl(r(sel), 25), pctl(r(sel), 75), median(x(sel))};
            else
                v = {0, NaN, NaN, NaN, NaN};
            end
            rows(end + 1, :) = [{R.blkSeasons{si}, R.blkDrivers{di}, full{j}}, v];          %#ok<AGROW>
        end
    end
end
writeCsv(fullfile(T, 'tab03_driver_anomalies_by_region.csv'), {'season', 'driver', 'region', 'cells', 'median_rho', ...
    'p25_rho', 'p75_rho', 'median_xi'}, rows);
nS = numel(R.blkSeasons);  nD = numel(R.blkDrivers);
if numel(G.keep) <= E.maxCsvCells
    [rC, cC] = ind2sub([G.HC G.WC], G.keep);
    [~, o] = sort((rC - 1) * G.WC + cC);                               % row-major order
    anyr = squeeze(any(any(isfinite(R.blkRho), 1), 2));
    o = o(anyr(o));
    M = [rC(o), cC(o), (G.origin(1) + (cC(o) - 0.5) * G.agg * G.cell) / 1e3, ...
         (G.origin(2) - (rC(o) - 0.5) * G.agg * G.cell) / 1e3, G.regC(o)];
    hdr = {'cell_row', 'cell_col', 'x_km', 'y_km', 'region'};
    fmt = '%d,%d,%.6g,%.6g,%d';
    for si = 1:nS
        for di = 1:nD
            M = [M, squeeze(R.blkRho(si, di, o)), squeeze(R.blkXi(si, di, o)), squeeze(R.blkN(si, di, o))];   %#ok<AGROW>
            sd = [R.blkSeasons{si} '_' R.blkDrivers{di}];
            hdr = [hdr, {['rho_' sd], ['xi_' sd], ['n_' sd]}];                              %#ok<AGROW>
            fmt = [fmt, ',%.6g,%.6g,%d'];                                                   %#ok<AGROW>
        end
    end
    fid = fopen(fullfile(T, 'tab03b_driver_anomalies_by_cell.csv'), 'w');
    fprintf(fid, '%s\n', strjoin(hdr, ','));
    fprintf(fid, [fmt '\n'], M');
    fclose(fid);
else                                                                    % too many cells for a CSV: maps
    nB = nS * nD;  Mr = nan(G.HC, G.WC, nB);  Mx = Mr;  Mn = Mr;  rows = cell(nB, 3);
    for si = 1:nS
        for di = 1:nD
            q = (si - 1) * nD + di;
            Mr(:, :, q) = toGridC(squeeze(R.blkRho(si, di, :)), G);
            Mx(:, :, q) = toGridC(squeeze(R.blkXi(si, di, :)), G);
            Mn(:, :, q) = toGridC(squeeze(R.blkN(si, di, :)), G);
            rows(q, :) = {q, R.blkSeasons{si}, R.blkDrivers{di}};
        end
    end
    writeCellTif(fullfile(T, 'tab03b_driver_rho_by_cell.tif'), Mr, G);
    writeCellTif(fullfile(T, 'tab03b_driver_xi_by_cell.tif'), Mx, G);
    writeCellTif(fullfile(T, 'tab03b_driver_pairs_by_cell.tif'), Mn, G);
    writeCsv(fullfile(T, 'tab03b_bands.csv'), {'band', 'season', 'driver'}, rows);
    clear Mr Mx Mn
end
% tab03c: driver-anomaly correlations by permafrost extent class of the cell
rows = {};
for si = 1:nS
    for di = 1:nD
        for k = 1:5
            rows(end + 1, :) = {R.blkSeasons{si}, R.blkDrivers{di}, k - 1, R.pfZones{k}, R.pfZoneN(si, di, k), ...
                R.pfZoneRho(si, di, k, 2), R.pfZoneRho(si, di, k, 1), R.pfZoneRho(si, di, k, 3), R.pfZoneXi(si, di, k)};   %#ok<AGROW>
        end
    end
end
writeCsv(fullfile(T, 'tab03c_driver_anomalies_by_permafrost.csv'), {'season', 'driver', 'permafrost_class', ...
    'permafrost_extent', 'cells', 'median_rho', 'p25_rho', 'p75_rho', 'median_xi'}, rows);
% tab04: rho and xi; tab05: Spearman matrices over pixel-days
rows = {};
for s = 1:3
    for j = 2:numel(R.vars)
        unit = 'pixel-days';  if j > 6, unit = 'pixels (season mean Albedo)'; end
        rows(end + 1, :) = {E.seasons{s, 1}, R.vars{j}, unit, R.rho{s}(j), R.xi{s}(j), R.xiN{s}(j)};   %#ok<AGROW>
    end
end
writeCsv(fullfile(T, 'tab04_dependence_rho_xi.csv'), {'season', 'covariate', 'unit', 'spearman_rho', 'chatterjee_xi', 'n'}, rows);
for s = 1:3
    writeCsv(fullfile(T, sprintf('tab05_spearman_pixeldays_%s.csv', E.seasons{s, 1})), [{''}, R.vars], ...
        [R.vars(:), num2cell(R.corr{s})]);
end
% tab09: five-year Feb-Oct mean Albedo by region (cells with every period) and its difference from the normal
rows = {};
for j = 1:10
    for q = 1:size(R.periods, 1)
        rows(end + 1, :) = {full{j}, sprintf('%d-%d', R.periods(q, 1), R.periods(q, 2)), R.perYears(q), R.perRegN(j), ...
            R.perReg(j, q), R.perRegBase(j), R.perReg(j, q) - R.perRegBase(j)};                           %#ok<AGROW>
    end
end
writeCsv(fullfile(T, 'tab09_five_year_means.csv'), {'region', 'period', 'years_read', 'cells', 'mean_Albedo_FebOct', ...
    'baseline_Albedo_FebOct', 'difference'}, rows);
end

function writeExtraTables(X, E, G)
% tables of stage 2b: tab04b, tab11, tab11b, tab12 and tab13 (every run)
full = [G.regNames(G.order), {'All Alaska'}];  T = E.tabDir;  yb = X.yb;  nMet = numel(X.metName);
% tab11: annual Albedo metrics of every region and year; tab11b: their trends
rows = cell(10 * numel(yb), 2 + nMet);  q = 0;
for j = 1:10
    for i = 1:numel(yb)
        q = q + 1;  rows(q, :) = [{full{j}, yb(i)}, num2cell(reshape(X.met(j, i, :), 1, []))];
    end
end
hdr = [{'region', 'year'}, regexprep(strcat(X.metName, '_', X.metUnit), '[^A-Za-z0-9]+', '_')];
writeCsv(fullfile(T, 'tab11_Albedo_annual_metrics.csv'), hdr, rows);
rows = {};
for j = 1:10
    for k = 1:nMet
        rows(end + 1, :) = {full{j}, X.metName{k}, X.metUnit{k}, X.metTrend.n(j, k), X.metTrend.slope(j, k), ...
            X.metTrend.tau(j, k), X.metTrend.p(j, k), X.metTrend.q(j, k)};                              %#ok<AGROW>
    end
end
writeCsv(fullfile(T, 'tab11b_Albedo_metric_trends.csv'), {'region', 'metric', 'unit', 'years', 'theil_sen_per_decade', ...
    'kendall_tau', 'mann_kendall_p', 'bh_fdr_q'}, rows);
% tab04b: rho and xi of Albedo vs the covariates by region
vars = {'NDSI', 'NDVI', 't2m', 'SWE', 'snowfall', 'elevation', 'slope', 'northness', 'latitude', 'peak NDVI', 'permafrost'};
rows = {};
for s = 1:3
    for j = 1:10
        for q = 1:numel(vars)
            unit = 'pixel-days';  if q > 5, unit = 'pixels (season mean Albedo)'; end
            rows(end + 1, :) = {E.seasons{s, 1}, full{j}, vars{q}, unit, X.dep.n(j, q, s), X.dep.rho(j, q, s), ...
                X.dep.xi(j, q, s)};                                                                     %#ok<AGROW>
        end
    end
end
writeCsv(fullfile(T, 'tab04b_dependence_rho_xi_by_region.csv'), {'season', 'region', 'covariate', 'unit', 'n', ...
    'spearman_rho', 'chatterjee_xi'}, rows);
% tab12: mean annual cycle of every driver by region (365-day calendar, 7-day smoothed)
rows = cell(10 * 365 * numel(X.drv), 5);  q = 0;
for v = 1:numel(X.drv)
    for j = 1:10
        for d = 1:365
            q = q + 1;  rows(q, :) = {X.drv{v}, X.drvUnit{v}, full{j}, d, X.drvCycMean(j, d, v)};
        end
    end
end
writeCsv(fullfile(T, 'tab12_driver_mean_annual_cycle.csv'), {'driver', 'unit', 'region', 'day_365', 'mean'}, rows);
% tab13: statewide daily Albedo and NDSI snow cover of the snow years
for k = 1:numel(X.snow)
    sw = X.snow(k);
    writeCsv(fullfile(T, sprintf('tab13_snow_year_%d_%d.csv', sw.y0, sw.y0 + 1)), {'date', 'Albedo', 'ndsi_snow_cover', ...
        'Albedo_7day', 'ndsi_7day', 'Albedo_baseline_mean', 'ndsi_baseline_mean', 'Albedo_valid_share', 'ndsi_valid_share'}, ...
        [cellstr(datestr(sw.dn, 'yyyy-mm-dd')), num2cell([sw.alb, sw.ndsi, sw.albS, sw.ndsiS, sw.albClim, sw.ndsiClim, ...
        sw.albCover, sw.ndsiCover])]);
end
fprintf('  wrote tab04b, tab11, tab11b, tab12 and %d tab13 table(s) to %s\n', numel(X.snow), T);
end

function writeCellTif(fp, A, G)
% HC x WC (x bands) map on the cell grid -> single-precision GeoTIFF (NaN = no value), with the
% grid and CRS of the Prepared static_ecoregion.tif of the same cells
tmpl = fullfile(G.gridDir, 'static_ecoregion.tif');
gi = georasterinfo(tmpl);
assert(isequal(gi.RasterSize(1:2), [G.HC G.WC]), '%s is not on the cell grid', tmpl);
gt = geotiffinfo(tmpl);
keys = gt.GeoTIFFTags.GeoKeyDirectoryTag;
if isfield(keys, 'Unknown'), keys = rmfield(keys, 'Unknown'); end     % geotiffwrite refuses unnamed keys (as alb_prepare.m)
assert(isfield(keys, 'ProjectedCSTypeGeoKey') && keys.ProjectedCSTypeGeoKey == 3338, '%s is not tagged EPSG:3338', tmpl);
geotiffwrite(fp, single(A), gi.RasterReference, 'GeoKeyDirectoryTag', keys);
end

function writeCsv(fp, hdr, rows)
% header (cell of names) and rows (cell matrix of numbers or text) -> CSV; NaN -> empty
fid = fopen(fp, 'w');  assert(fid > 0, 'Cannot write %s', fp);
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(rows, 1)
    c = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if ischar(v)
            if any(v == ','), v = ['"' v '"']; end
            c{j} = v;
        elseif ~isfinite(v)
            c{j} = '';
        elseif v == round(v) && abs(v) < 1e15
            c{j} = sprintf('%d', round(double(v)));             % counts in full
        else
            c{j} = sprintf('%.6g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(c, ','));
end
fclose(fid);
end

% ----------------------------- file access (as alb_prepare.m) -----------------------------

function Fk = openTiff(fp)
t = Tiff(fp, 'r');
info = imfinfo(fp);
Fk = struct('t', t, 'fp', fp, 'nB', info(1).SamplesPerPixel, 'H', info(1).Height, 'W', info(1).Width, ...
            'tw', getTag(t, 'TileWidth'), 'th', getTag(t, 'TileLength'));
end

function closeTiff(Fk)
close(Fk.t);
end

function a = readTile(Fk, rr, cc)
% One tile, all bands -> pixels (column-major) x bands. Falls back to imread if needed.
try
    a = readEncodedTile(Fk.t, computeTile(Fk.t, [rr(1) cc(1)]));
    a = a(1:numel(rr), 1:numel(cc), :);
catch
    a = imread(Fk.fp, 'PixelRegion', {[rr(1) rr(end)], [cc(1) cc(end)]});
end
assert(size(a, 3) == Fk.nB, 'Band count mismatch reading %s', Fk.fp);
a = reshape(a, [], size(a, 3));
end

function dn = tagDates(fp)
% band dates (datenum) from the GDAL band descriptions 'dYYYYMMDD'
tok = regexp(gdalMetadata(fp), 'sample="(\d+)"[^>]*role="description">d(\d{8})<', 'tokens');
assert(~isempty(tok), '%s: no band dates (GDAL band descriptions dYYYYMMDD) found', fp);
smp = cellfun(@(c) str2double(c{1}), tok);  v = cellfun(@(c) str2double(c{2}), tok);
[smp, k] = sort(smp);  v = v(k);
assert(isequal(smp, 0:numel(smp) - 1), '%s: band descriptions do not cover bands 1..%d', fp, numel(smp));
dn = datenum(floor(v / 1e4), mod(floor(v / 100), 100), mod(v, 100));
dn = dn(:);
assert(all(diff(dn) > 0), '%s: band dates are not strictly increasing', fp);
info = imfinfo(fp);
assert(numel(dn) == info(1).SamplesPerPixel, '%s: %d band dates for %d bands', fp, numel(dn), info(1).SamplesPerPixel);
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
