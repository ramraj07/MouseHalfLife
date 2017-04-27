function mouseHalfLife
%This function is effectively a script which is used to fit pharmacokinetic
%data to various exponential decay-models in order to find the half-life of
%various biomolecules in mice. Input is the data from all the counts (or
%ELISA values) and must be specified as a CSV (Comma Seperated Value) file
%in the prescribed format. Output is a PDF with all the fitting data and
%calculated half-lives.

dataFilename= 'lsfdata.csv';
%Enter the name of the csv file that contains the half-life data in the
%prescribed format. CSV is a special data-only format you can save in excel; the
%file should contain only a single table of the form:
%
%     Group:    loxp+/+     loxp+/+     loxp+/+  ...
%     Mouse:    138M-       138Mr       138Ml       ...
%     0         75.4        78.6        94.3        ...
%     8         64.8        66          83.4        ...
%     10        61.6        63          78.6        ...
%     14        58.4        58.2        71.5        ...
%     24        49.3        48.7        58.6        ...
%     34        42.6        39.8        45.5        ...
%     54        30          27.2        30.3        ...
%     71        24.4        20.6        24.4        ...
%     89        18.68       15.18       17.27       ...
%     110       13.44       11.08       11.76       ...
%     ...       ...         ...         ...         ...
%     ...       ...         ...         ...         ...
%     ...       ...         ...         ...         ...
%
% The first column contains the timepoints at which readings were taken and
% the first two rows tell the group and mouse identification information.
% Give exactly same text for mice which are in the same group, thats how
% they are "grouped" for the report.


ExperimentTitle = 'Pka in KO and GKO 03.27.10';
ExperimentDescription = '';
% Provide description of the experiment, be as elaborate as you want

ExperimentDate = '03-27-10';
ExperimenterName=  '';

CorrectBackground = 'no';
% If there is background correction to be done, enter 'no'. Otherwise enter
% 'yes'.

IsotopeHalfLife =1430.4;
% In hours (or whatever time unit used in the data

TimeUnits = 'hours';

AutomaticInitialConditions = 'yes';
% If you want initial conditions to be determined automatically, put 'yes'.
% Otherwise put 'no' and provide values for the conditions below. If
% single-exponential fitting method is used only the HalfLife1 and Coeff1
% values will be utilized.

InitialConditions_HalfLife1 = 20;
InitialConditions_HalfLife2 = 220;
InitialConditions_Coeff1 = 0.03;
InitialConditions_Coeff2 = 0.003;


% BEST DEFAULT FITTING METHOD: 2
fittingMethod = 1;
% There are three fitting methods you can use, and you choose by putting
% 1,2 or 3 above. The explanations for them are:
%
% 1 = Independent-coefficient biexponential fit
%       In this method the equation used is,
%           y = C1*exp(-k1*t)+C2*exp(-k2*t)
%       This is the traditionally used method, where you expect the curve
%       to follow a biexponential curve. k1 and k2 are the half lives.
%
% 2 = Semi-dependent-coefficient biexponential fit
%       The equation used in this method is,
%           y = 0.001*(C1*C2*exp(-k1*t)+C1*(100-C2)exp(-k2*t))
%       While this equation might look different, it means exactly the same
%       as the Independent-coefficient equation; the only difference is
%       that here the Constants C1 and C2 indicate physiologically relevant
%       values, viz. C1 is the approximate starting value at which the
%       readings start and C2 is the % contribution of the first
%       exponential to the whole decay process. The results obtained from
%       (1) or (2) are generally similar upto the third digit (after that
%       the decimals change a bit because of nuances in the fitting
%       landscape)
%
% 3 = Single exponential fit
%       The equation used here is,
%           y = C*exp(-k*t)
%       This should be used in special cases where the decay looks more
%       like a single exponential one rather than biexponential

% % % % % % DependentExponentials = 'yes';

logscalefitting = true;

% Fitting Paramters (generally don't have to meddle with these at all)
MaxFunEvals = 100000;
MaxIter = 10000;
TolFun = 10e-8;
TolX = 10E-9;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% END OF USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
datestrval = datestr(now,'mmmm-dd-yyyy.HH.MM.SS');
%% Reading and parsing of the data
%read the data from the csv file
data = importdata(dataFilename,',',2);
timePoints = data.data(:,1);
noOfTimePoints = length(timePoints);
miceNames = data.colheaders(2:end);
noOfMice = length(miceNames);
groupNames = cell(noOfMice,1);

% isolate and find unique group names
gtext = [' ',data.textdata{1}];
[temp gtext] = strtok(gtext,',');
for i=1:noOfMice
    [temp gtext] = strtok(gtext,',');
    groupNames{i} = temp;
end
readings = data.data(:,2:end);

% Find unique groups and groupId of each mouse
clear groups
groupId = zeros(noOfMice,1);
noOfGroups=1;
groups{noOfGroups} = groupNames{1};
groupId(find(strcmp(groupNames,groups{noOfGroups}))) = noOfGroups;
for i=2:length(groupNames)
    if groupId(i)~=0
        continue;
    end
    noOfGroups=noOfGroups+1;
    groups{noOfGroups} = groupNames{i};
    groupId(find(strcmp(groupNames,groups{noOfGroups}))) = noOfGroups;
end
if ~strcmpi('mouse:',(data.colheaders{1}))
    error(['Second column is not labeled ''Mouse:''! Please adhere to',...
        'format shown above and READ ALL THE INSTRUCTIONS!']);
end
if ~strcmpi('group:',data.textdata{1}(1:6))
    error(['First column is not labeled ''Group:''! Please adhere to',...
        'format shown above and READ ALL THE INSTRUCTIONS!']);
end
%% Correct Data for various effects

% Correct for background
if strcmp(CorrectBackground,'no')~=1
    error('Background Correction option has not yet been incorporated in this version of the program.');
end
bgCorrectedData = readings;



% Correct for isotope decay
% The correction is done with the formula A0 = e^(k0*t)*At where A0 is what
% the corrected count would be, k0 is the decay constant (-ln2 / half life)
% At is the measurement done at time t.

k0 = log(0.5) / (-IsotopeHalfLife);
decayCorrectedData = zeros(size(readings,1),size(readings,2));
bgCorrectedNormalizedData = zeros(size(readings,1),size(readings,2));
meanNormalized = zeros(noOfTimePoints,noOfGroups);
for i=1:noOfMice
    %normalize data for exporting and plotting
    bgCorrectedNormalizedData(:,i) = 100*bgCorrectedData(:,i)/max(bgCorrectedData(:,i));
end
for i=1:noOfTimePoints
    decayCorrectedData(i,:) = exp(k0*timePoints(i))*bgCorrectedData(i,:);
    
end
for i=1:noOfGroups
    for j=1:noOfTimePoints
        valuess = bgCorrectedNormalizedData(j,groupId==i);
        meanNormalized(j,i) = mean(valuess(~isnan(valuess)));
    end
end

ans=questdlg('You are attempting to run the fit with highly experimental and non-recommended log-space curve fitting option. Are you aware of its implications and want to proceed?',...
    'Highly Expereimental Method','No I''m not aware, don''t proceed',...
    'Yes I know what I''m doing, proceed','No I''m not aware, don''t proceed');
if ~strcmp(ans,'Yes I know what I''m doing, proceed')
    error('asfss');
end


%% Determine Initial Conditions if required or parse them

InitialConditions = cell(noOfMice,1);
switch(AutomaticInitialConditions)
    case 'yes'
        for i=1:noOfMice
            % ymax is the maximum initial starting point of the decay curve
            dataForThisMouse = decayCorrectedData(:,i);
            ymax =  max(decayCorrectedData(:,i));
            % xmid is the time at which approximately y reaches ymax/2
            %
            [ onebefore] = find(dataForThisMouse<(ymax/2),1,'first');
            if onebefore == length(dataForThisMouse)
                onebefore = onebefore-1;
            end
            xmid=polyval(...
                polyfit(...
                timePoints(onebefore:onebefore+1),...
                dataForThisMouse(onebefore:onebefore+1),...
                1)...
                ,ymax/2);
            %         xmid = ((max(timePoints)-min(timePoints))/2);
            
            
            switch(fittingMethod)
                case 1
                    % Independent exponentials; the first exp is assumed to
                    % contribute to 2/3rd of the effect as a starting
                    % position, so C1 starts with 2/3 of ymax and the rate
                    % const is assumed to be around 0.3 / the timepoint at
                    % which half the decay has occurred and the other
                    % exponential vice versa
                    %
                    InitialConditions{i} = [...
                        ymax*(2/3); ...
                        0.3 / xmid;...
                        ymax/3;...
                        3 / xmid ];
                case 2
                    % Similar method as above but the C1 is taken to be
                    % near ymax and C2 is taken to be 95 (which means 95%
                    % of the decay is due to the first exponential, a case
                    % reasonable to most WT mice IgG decay curves
                    InitialConditions{i} = [...
                        ymax; ...
                        0.3 / xmid;...
                        95;...
                        3 / xmid ];
                case 3
                    InitialConditions{i} = [...
                        ymax;...
                        0.2 / xmid];
            end
        end
    case 'no'
        for i=1:noOfMice
            switch(fittingMethod)
                case  {1,2}
                    InitialConditions{i} = [ log(2)/InitialConditions_HalfLife1, ...
                        InitialConditions_Coeff1 ,...
                        log(2)/InitialConditions_HalfLife2,...
                        InitialConditions_Coeff2 ];
                case 3
                    InitialConditions{i} = [ log(2)/InitialConditions_HalfLife1, ...
                        InitialConditions_Coeff1];
            end
            
        end
end

%% Run the fitting for each mouse
OptOptions = ...
    optimset('LargeScale','off','LevenbergMarquardt','on','Diagnostics','on', ...
    'MaxFunEvals', MaxFunEvals,'MaxIter',MaxIter,'TolFun',TolFun,'TolX',TolX,'Display','iter') ;
weights = ones(length(timePoints),1);
results = zeros(noOfMice,4);
HalfLives = zeros(noOfMice,2);
simulations = decayCorrectedData*0;
simy1 = simulations;
simy2 = simulations;

maxVal = max(max(decayCorrectedData));
maxk = log(2)/1;
mink = log(2)/max(timePoints);
resultsCell = cell(noOfMice,5);

%lower Bound / Upper Bound specification
switch(fittingMethod)
    case 1
        lb = [0 1e-5 0 .01];
        ub = [maxVal 0.5 maxVal 0.5];
    case 2
        lb = [1 log(2)/1000 0 log(2)/1000];
        ub = [maxVal*2 0.5 100 0.5];
    case 3
        lb = [1 1e-5];
        ub = [maxVal*2 0.5];
end

disp('Fitting the data...');
for i=1:noOfMice
    
    
    
    options = optimset('maxiter', MaxIter,...
        'levenbergmarquardt', 'off',...
        'maxfunevals', MaxFunEvals,...
        'display', 'off',...
        'tolx', TolX,...
        'tolfun', TolFun);
    
    %check for NaNs and remove those points
    if isempty(find(isnan(decayCorrectedData(:,i))))
        prunedTimepoints = timePoints;
        dataGoingIn = decayCorrectedData(:,i);
    else
        nonNaNpoints = ~isnan(decayCorrectedData(:,i));
        prunedTimepoints = timePoints(nonNaNpoints);
        dataGoingIn = decayCorrectedData(nonNaNpoints,i);
    end
    if logscalefitting
        
        switch(fittingMethod)
            case 1
                [parameters,resnorm]=lsqcurvefit(@decay_biexp_log,...
                    InitialConditions{i},prunedTimepoints,log(dataGoingIn),...
                    lb,ub,options);
            case 2
                [parameters,resnorm]=lsqcurvefit(@decay_biexp_dep_log,...
                    InitialConditions{i},prunedTimepoints,log(dataGoingIn),...
                    lb,ub,options);
            case 3
                [parameters,resnorm]=lsqcurvefit(@decay_exp_log,...
                    InitialConditions{i},prunedTimepoints,log(dataGoingIn),...
                    lb,ub,options);
        end
        
    else
        switch(fittingMethod)
            case 1
                [parameters,resnorm]=lsqcurvefit(@decay_biexp,...
                    InitialConditions{i},prunedTimepoints,dataGoingIn,...
                    lb,ub,options);
            case 2
                [parameters,resnorm]=lsqcurvefit(@decay_biexp_dep,...
                    InitialConditions{i},prunedTimepoints,dataGoingIn,...
                    lb,ub,options);
            case 3
                [parameters,resnorm]=lsqcurvefit(@decay_exp,...
                    InitialConditions{i},prunedTimepoints,dataGoingIn,...
                    lb,ub,options);
        end
    end
    
    if length(parameters)==2
        parameters = [parameters' NaN NaN]';
    end
    results(i,:) = parameters;
    
    
    %extract parameters
    
    coeff1 = parameters(1);
    kd1    = parameters(2);
    coeff2 = parameters(3);
    kd2    = parameters(4);
    
    %calculate simulated fits
    switch(fittingMethod)
        case 1
            simy1(:,i) = coeff1*exp(-kd1*timePoints);
            simy2(:,i) = coeff2*exp(-kd2*timePoints);
        case 2
            simy1(:,i) = coeff1*coeff2*0.01*exp(-kd1*timePoints);
            simy2(:,i) = coeff1*(100-coeff2)*0.01*exp(-kd2*timePoints);
        case 3
            simy1(:,i) = coeff1*exp(-kd1*timePoints);
    end
    
    simulations(:,i) =simy1(:,i)+simy2(:,i);
    halflives(i,1) = log(2)/ parameters(2);
    halflives(i,2) = log(2)/ parameters(4);
    resultsCell{i,1}= i;
    resultsCell{i,2}= groups{groupId(i)};
    resultsCell{i,3}= miceNames{i};
    resultsCell{i,4}=  halflives(i,1);
    resultsCell{i,5}=  halflives(i,2);
end
residuals = decayCorrectedData-simulations;
stdres = std(residuals,0,1);

fprintf(1,'Generating report...');
%% Generate the report
close all

maxchars = 50;

% s = strtrim(ExperimentDescription);
% tokens = regexp(s, sprintf('(\\S\\S{%d,}|.{1,%d})(?:\\s+|$)',...
%     maxchars, maxchars), 'tokens').';
% get_contents = @(f) f{1};
% c = cellfun(get_contents, tokens, 'UniformOutput', false);
% c = deblank(c);
switch(fittingMethod)
    case 1
        FitEquation = 'y = C_1{\ite}^{-\itk_1t}+C_2{\ite}^{-\itk_2t}';
    case 2
        FitEquation = 'y = {0.01}({C_1}{C_2}{\ite}^{-\itk_1t}+{C_1}({100-C_2}){\ite}^{-\itk_2t})';
    case 3
        FitEquation = 'y = C{\ite}^{-\itkt}';
end

figuretext = [{'\fontsize{16}\bfMouse Half Life analysis\rm\fontsize{14}'},...
    linewrap([ExperimentTitle ' '],25)',...
    {['\rm\fontsize{12}Experimenter Name: ',ExperimenterName],...
    ['Performed on ',ExperimentDate],...
    ' ',...
    ['Isotope half-life: ',num2str(IsotopeHalfLife),' ',TimeUnits],...
    ['Determine Initial Conditions automatically? \bf',...
    AutomaticInitialConditions],...
    'Fit equation:\rm',...
    ' ',...
    FitEquation,...
    ' ',...
    'Groups in the study\rm',...
    strvcat(groups{:}),...
    ' ',...
    '\bfExperiment description\rm\fontsize{10}'},...
    linewrap([ExperimentDescription,' '],30)'];

% for i=1:length(c)
%     figuretext{length(figuretext)+1} = c{i};
% end
figuretext{end+1} = '';
figuretext{end+1} = ['Analysis done @ ',datestrval];
h1= figure('Position',[0 0 1280 896],'units','normalized');
if logscalefitting
    annotation('textbox',[0.025 .025 .9 .05],'string',...
        'Highly experimental Log-scale based fitting method used; data should not be used for further analysis!',...
        'fontsize',16,'color','r')
end
h2= subplot(1,3,3);
text (0.5,0.5,figuretext,'HorizontalAlignment','Center',...
    'VerticalAlignment','middle','fontsize',10);
set(gca,'visible','off');
%f = figure('Position',[200 200 960 678]);
dat =  resultsCell;
columnname =   {'', 'Group Name','Mouse', 'Half Life 1', 'Half Life 2'};
columnformat = {'numeric', 'char', 'numeric', 'numeric'};
columneditable =  [false false true true];
h3 = subplot(1,3,[1 2]); poss = get(h3,'position');delete(h3);
t = uitable( 'Data', dat,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'units','normalized',...
    'RowName',[],...
    'Position',poss,.....
    'columnwidth','auto',...
    'FontSize',9      );
set(gcf, 'PaperOrientation', 'landscape');
set(gcf,'PaperPosition', [0.5 0.5 10 7.5]);
%%
print(h1,'resultstemp','-dpsc');
%% Plot averaged plot and also the raw data
clf
colorCell = {'r','g','b','m','c','k','y',...
    'r','g','b','m','c','k','y',...
    'r','g','b','m','c','k','y'};
for i = 1:length(groups)
    plot(timePoints,...
        meanNormalized(:,i),...
        '-o',    'color',colorCell{i},'linewidth',2.5,...
        'MarkerFaceColor',colorCell{i},...
        'markersize',4);
    
    hold on;
end
for i = 1:length(groups)
    
    errorbar(timePoints,...
        meanNormalized(:,i),...
        std(bgCorrectedNormalizedData(:,groupId==i),[],2),...
        'linestyle','none','color',colorCell{i},...
        'linewidth',0.5);
    hold on
end
hold off
grid on
legend(groups);
axis([0 max(timePoints) 0 100]);
if logscalefitting
    annotation('textbox',[0.025 .025 .9 .05],'string',...
        'Highly experimental Log-scale based fitting method used; data should not be used for further analysis!',...
        'fontsize',16,'color','r')
end
xlabel(['Time (', TimeUnits,')']);
ylabel('Background corrected normalized values (% initial)');
title('\fontsize{14}\bfPlot of grouped, background subtracted normalized decay curves');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf,'PaperPosition', [0.5 0.5 10 7.5]);

print(h1,'resultstemp','-dpsc','-append');
%%
clf
%axes('xtick',[],'ytick',[])
origdata = {'{\bfOriginal Data analyzed:}','',...
    ['Timepoints: ',num2str(timePoints','%d,'),' ',TimeUnits]};

for i=1:noOfMice
    origdata{end+1} = ['{\bfMouse ', miceNames{i},'}: Group ',groupNames{i}];
    origdata{end+1} = ['Readings: ',num2str(readings(:,i)','%0.0f,')];
end
at = annotation('textbox',[0.05 .05 .9 .9],'string',...
    origdata,'fitboxtotext','on','linestyle','none');
if logscalefitting
    annotation('textbox',[0.025 .025 .9 .05],'string',...
        'Highly experimental Log-scale based fitting method used; data should not be used for further analysis!',...
        'fontsize',16,'color','r')
end
% poss = get(at,'position');
% if poss(1)+poss(3)>1 || poss(2)+poss(4)>1
%     %the box is going too far!
%     scaledown= max(poss(1)+poss(3)/1 , poss(2)+poss(4))
set(gcf, 'PaperOrientation', 'landscape');
set(gcf,'PaperPosition', [0.5 0.5 10 7.5]);
%%
print(h1,'resultstemp','-dpsc','-append');

%%
% plot each mouse
disp('   ');
try
    resultsFolder=['Results_',ExperimentTitle,'_',datestrval];
    mkdir(resultsFolder);
catch e
    resultsFolder=['Results_',datestrval];
    mkdir(resultsFolder);
end
fp = fopen([resultsFolder,filesep,'results.csv'],'w');
if isempty(ferror(fp))
    fprintf(fp,'Group,Mouse,Half Life 1 (%s), Half Life 2 (%s)',...
        TimeUnits,TimeUnits);
end
for i=1:noOfMice
    
    
    clf
   
    subplot(4,2,[1 3]);
    
    plot(timePoints,decayCorrectedData(:,i),'-ro','LineWidth',2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    ax = axis;
    title(['\fontsize{15}',groups{groupId(i)},' | Mouse: ',miceNames{i}]);
    ylabel('Adjusted Counts');
    
    subplot(4,2,[2 4]);
    plot(timePoints,simulations(:,i),'-','LineWidth',2      );
    hold on;
    plot(timePoints,decayCorrectedData(:,i),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    plot(timePoints,simy1(:,i),'--');
    plot(timePoints,simy2(:,i),':');
    axis(ax);
    ylabel('Adjusted Counts');
    title(['\fontsize{12}Results of NLR analysis']);
    legend('Simulation curve','Adjusted original data',...
        'Exponential 1 only','Exponential 2 only');
    grid(gca,'on');
    subplot(4,2,[5 7]);
    plot(timePoints,simulations(:,i),'-','LineWidth',2);
    set(gca,'yscale','log');
    axis(ax);
    hold on
    plot(timePoints,decayCorrectedData(:,i),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    legend('Simulated data','Adjusted original data');
    title('\fontsize{12}Simulation+data on semi-log');
    grid(gca,'on');
    subplot(4,2,6);
    plot(timePoints,abs(residuals(:,i)),'o');
    ylabel('Adjusted Counts'); xlabel(['Time (',TimeUnits,')']);
    title('\fontsize{12}Residuals');
    
    subplot(4,2,8);
    set(gca,'visible','off');
    
    text1{1} = ['\fontsize{12}Half Life 1: \bf',num2str(halflives(i,1),4),...
        ' ',TimeUnits];
    text1{2} = ['\rm\fontsize{12}Half Life 2: \bf',...
        num2str(halflives(i,2),4),' ',TimeUnits];
    text1{3} = '\rm\fontsize{12}Fit Parameters:\rm';
    
    text1{4} = ['\fontsize{10}Coefficient 1: ',num2str(results(i,1),4),...
        '; Coefficient 2: ',...
        num2str(results(i,3),4)];
    text1{5} = ['\fontsize{10}k1: ',num2str(results(i,2)),'; k2: ',...
        num2str(results(i,4));];
    text1{6} = '\fontsize{12}Initial Conditions:\fontsize{10}';
    if fittingMethod == 3
        text1{7} = sprintf('C1: %.2f; k1: %.4f;',...
            InitialConditions{i}(1),InitialConditions{i}(2));
    else
        text1{7} = sprintf('C1: %.2f; k1: %.4f; C2: %.2f; k2: %.4f;',...
            InitialConditions{i}(1),InitialConditions{i}(2),...
            InitialConditions{i}(3),InitialConditions{i}(4));
    end
    text(0,.5,text1);
     if logscalefitting
        annotation('textbox',[0.025 .025 .9 .05],'string',...
            'Highly experimental Log-scale based fitting method used; data should not be used for further analysis!',...
            'fontsize',16,'color','r')
    end
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf,'PaperPosition', [0.5 0.5 10 7.5]);
    fprintf(1,'%02d. %s\t%s\t%03.2f %s\t\t%03.2f %s\n',i,...
        groups{groupId(i)},miceNames{i},halflives(i,1),TimeUnits,...
        halflives(i,2),TimeUnits);
    if isempty(ferror(fp))
        fprintf(fp,'%s,%s,%f,%f\n',i,...
            groups{groupId(i)},miceNames{i},halflives(i,1),halflives(i,2));
    end
    %%
    print(gcf,'resultstemp','-dpsc','-append');
    
    
end
try
    fclose(fp);
end
close

% Make a "Results" folder which is unique for this execution by putting a
% timestamp also, and then copy this script file, save a mat file of all
% the variables' values and also a pdf of the results. Also copy the
% original csv file that was used.
copyfile([mfilename,'.m'],[resultsFolder,filesep,datestrval,'_script.m']);
copyfile(dataFilename,resultsFolder);
save([resultsFolder,filesep,datestrval,'_analysisdata.mat']);
pdfFile = [resultsFolder,filesep,datestrval,'_results.pdf'];
ps2pdf('psfile','resultstemp.ps','pdffile',pdfFile);
delete('resultstemp.ps');
open(pdfFile);
mkdir([resultsFolder,filesep,'generated CSV files'])
%%
fprintf(1,'Making CSV files...');

fp = fopen([resultsFolder,filesep,'generated CSV files',filesep,'NormalizedData.csv'],'w');
%fprintf(fp,'Group,Mouse,\n');
fprintf(fp,',Time (%s):,%s\n',TimeUnits,num2str(timePoints','%d,'));
for i=1:noOfMice
    numbersstring = strrep(num2str(bgCorrectedNormalizedData(:,i)','%f,'), 'NaN', '');
    
    fprintf(fp,'%s,%s,%s\n',groupNames{i},miceNames{i},numbersstring);
end
fprintf(fp,'\n\nGroup-wise averages\nTime (%s):,,%s\n',TimeUnits,num2str(timePoints','%d,'));
for i=1:noOfGroups
    numbersstring = strrep(num2str(...
        meanNormalized(:,i)','%f,'),...
        'NaN', '');
    
    fprintf(fp,'%s,,%s\n',groups{i},numbersstring);
    try
        fp2 = fopen([resultsFolder,filesep,'generated CSV files',filesep,...
            groups{i},'.csv'],'w');
        locs = find(groupId==i);
        for j=1:noOfTimePoints
            fprintf(fp2,'%d %s,',timePoints(j),TimeUnits);
        end
        for j=1:length(locs)
            numbersstring = strrep(num2str(decayCorrectedData(:,locs(j))','%f,'), 'NaN', '');
            fprintf(fp2,'%s\n',numbersstring);
        end
        fclose(fp2);
    catch e
        try
            if exists('fp2','var')==1
                if fp2==-1
                    warning(['Could not create csv file ',groups{i},'.csv']);
                end
                fclose(fp2);
            else
                warning(['Could not write to csv file ',groups{i},'.csv']);
            end
        end
    end
end

fclose(fp);



fprintf(1,'Done.');
end
% various fitting functions
function F = decay_biexp(x,xdata)
F = x(1)*exp(-x(2).*xdata)+ x(3)*exp(-x(4).*xdata);
end
function F = decay_biexp_dep(x,xdata)
F = x(1)*x(3)*0.01*exp(-x(2).*xdata)+ x(1)*(100-x(3))*0.01*exp(-x(4).*xdata);
end
function F = decay_exp(x,xdata)
F = x(1)*exp(-x(2).*xdata);
end
% various fitting functions
function F = decay_biexp_log(x,xdata)
F = log(x(1)*exp(-x(2).*xdata)+ x(3)*exp(-x(4).*xdata));
end
function F = decay_biexp_dep_log(x,xdata)
F = log(x(1)*x(3)*0.01*exp(-x(2).*xdata)+ x(1)*(100-x(3))*0.01*exp(-x(4).*xdata));
end
function F = decay_exp_log(x,xdata)
F = log(x(1)*exp(-x(2).*xdata));
end

%%

%This is ps2pdf, a method used by the program to generate the PDF file from
%the postscript files made. This requires ghostscript be installed in the
%machine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ps2pdf(varargin)
%PS2PDF Function to convert a PostScript file to PDF using Ghostscript
%
%  Converts a postscript file into PDF. The resulting PDF file will contain
%  one page for each page defined in the postscript files, so a multi-page
%  postscript file, like those generated by using the '-append' option of
%  MATLAB's print command, can be used to generate a multi-page PDF file.
%
%   Ghostscript is a third-party application currently supplied with
%   MATLAB. The caller may also specify a different version of Ghostscript
%   to use.
%
%   PS2PDF expects to be called with a set of parameter-value pairs. The
%   order of these is unimportant, but required parameters MUST be
%   specified
%
%   In the list below, required parameters are marked with an asterisk *
%        NOTE: ps2pdf can not use MATLAB's version of Ghostscript
%              in a deployed application; you MUST provide a
%              the path to a separate instance of Ghostscript. This
%              parameter is marked with a double asterisk **
%
%   Parameter:        Value:
%   *  psfile         full or relative path to the postscript file to convert
%
%   *  pdffile        full or relative path to the pdf file to create
%
%   ** gscommand      path to Ghostscript executable to use; this will try
%                     to default to the version of Ghostscript shipped with
%                     MATLAB, if any. If this value is specified you should
%                     also specify the gsfontpath and gslibpath values.
%
%                     ** See note on deployed applications, above.
%
%      gsfontpath     full path to the Ghostscript font files
%                     If a gscommand is specified then this path should
%                     also be specified and reference the same Ghostscript
%                     version
%
%      gslibpath      full path to the Ghostscript library (.ps) files.
%                     If a gscommand is specified then this path should
%                     also be specified and reference the same Ghostscript
%                     version
%
%                     If gscommand is NOT specified and we can determine
%                     the version of Ghostscript, if any, shipped with
%                     MATLAB, then this value will be overridden to use the
%                     path that references MATLAB's version of Ghostscript
%
%      gspapersize    paper size to use in the created .pdf file. If not
%                     specified or the specified value is not recognized
%                     it will use whatever default paper size is
%                     built into the version of Ghostscript being run
%
%                         NOTE: no scaling of the input occurs - it's simply
%                         placed on a page with the specified paper size.
%
%                         Valid values for gspapersize are:
%                              'letter', 'ledger', 'legal', '11x17',
%                              'archA', 'archB', 'archC', 'archD', 'archE',
%                              'a0', 'a1', 'a2', 'a3','a4', 'a5',
%                              'a6', 'a7', 'a8', 'a9', 'a10'
%
%      deletepsfile   0 to keep the input ps file after creating pdf
%                     non-zero to delete the input ps file after creating pdf
%                     Default is 0: keep the input ps file (do NOT delete it)
%                        NOTE: if the pdf creation process fails, the input
%                        PS file will be kept regardless of this setting
%
%      verbose        0 to suppress display of status/progress info;
%                     non-zero to allow display of status/progress info
%                     Default is 0 (no display)
%
% Example usage:
%    use MATLAB's version of Ghostscript to generate an A4 pdf file
%      ps2pdf('psfile', 'input.ps', 'pdffile', 'output.pdf', 'gspapersize', 'a4')
%
%    use a local copy of Ghostcript to generate a file, and display some
%    status/progress info while doing so.
%      ps2pdf('psfile', '../reports/input.ps', 'pdffile', 'c:\temp\output3.pdf', ...
%            'gspapersize', 'a4', 'verbose', 1, ...
%            'gscommand', 'C:\Program Files\GhostScript\bin\gswin32c.exe', ...
%            'gsfontpath', 'C:\Program Files\GhostScript\fonts', ...
%            'gslibpath', 'C:\Program Files\GhostScript\lib')
%
%    use MATLAB's version of Ghostscript to generate a pdf file and delete
%    the input.ps file when done
%      ps2pdf('psfile', 'input.ps', 'pdffile', 'output.pdf', 'gspapersize', 'a4', 'deletepsfile', 1)

%   Update log:
%      Jun 16, 2010: added check for deployed application
%      May 19, 2010: wrapped filenames sent to Ghostscript in quotes
%      May 06, 2010: updated how Ghostscript is found, don't rely on MATLAB version #
%      Aug 15, 2008: fixed bug where embedded space in location of
%                    MATLAB's version of Ghostscript caused ps2pdf to fail.
%      Apr 16, 2008: added deletepsfile option

%   Copyright 2008-2010 The MathWorks, Inc.

if nargin < 1
    error('ps2pdf:parameters', 'No parameters specified. Type ''help ps2pdf'' for details on how to use this function.');
end

% parse input args
gsData = LocalParseArgs(varargin{:});

% setup the file that tells GS what we want it to do
gsData = LocalCreateResponseFile(gsData);

gsDebug = 0;
if gsData.verbose
    fprintf('ps2pdf: input settings are:\n');
    if isfield(gsData, 'paperSizes')
        gsData = rmfield(gsData, 'paperSizes');
    end
    gsData  %#ok<NOPRT>
    fprintf('ps2pdf: response file for Ghostscript is:\n');
    type(gsData.responseFile);
    gsDebug = 1;
end

%to hold results/status from system call
s = 0; %#ok<NASGU>
r = ''; %#ok<NASGU>

% run Ghostscript to convert the file
if gsData.useBuiltin
    [s, r] = gsData.cmd(['@' gsData.responseFile], gsData.psFile, gsDebug);
else
    [s, r] = system([gsData.cmd ' @"' gsData.responseFile '" "' gsData.psFile '"']);
end

if gsData.verbose
    disp( ['Ghostscript STDOUT: ' num2str(s) ] );
    disp( ['Ghostscript STDERR: ' r ] );
else
    delete(gsData.responseFile)
end

if s && ~isempty(r)
    error('ps2pdf:ghostscript',  ['Problem converting PostScript. System returned error: ' num2str(s) '.' r])
elseif s
    error('ps2pdf:ghostscript',  ['Problem calling GhostScript. System returned error: ' num2str(s)])
end

%if after all this we still couldn't create the file, report the error
fid = fopen( gsData.pdfFile, 'r');
if ( fid == -1 )
    error('ps2pdf:ghostscript', '%s', [ 'Ghostscript could not create ''' gsData.pfdFile '''.' ])
else
    fclose( fid );
end

% if we get here, we successfully created pdf file; delete ps file if
% requested to do so
if gsData.deletePSFile
    delete(gsData.psFile);
end

end

%local function to parse arguments and fill in the gsData structure
%  .psFile - postscript file to convert
%  .pdfFile - pdf file to create
%  .cmd     - path/name of Ghostscript command to run or handle to gscript
%             builtin
%  .useBuiltin - 1 if using builtin gs command
%  .fontPath - path to the Ghostscript fonts, if any
%  .libPath - path to the Ghostscript libs (Ghostscript .ps files), if any)
%  .paperSize - paper size to set for resulting .pdf file
%  .deletePSFile - 0 to keep (not delete) the input ps file if pdf created ok
%  .verbose - if non-zero, display some status/progress info to command window
function gsData = LocalParseArgs(varargin)
gsData.paperSizes = {'letter', 'ledger', 'legal', '11x17', 'archA', 'archB', ...
    'archC', 'archD', 'archE', 'a0', 'a1', 'a2', 'a3','a4', 'a5', ...
    'a6', 'a7', 'a8', 'a9', 'a10'};

%default values for some settings
gsData.verbose      = 0;
gsData.useBuiltin   = 0;
gsData.deletePSFile = 0;

for i = 1 : 2 : length(varargin)-1
    param_arg = varargin{i};
    param_value = varargin{i+1};
    switch(lower(param_arg))
        % path to ps file to conver
        case 'psfile'
            if ~exist(param_value, 'file')
                error('print:ghostscript', ...
                    'Can not find postscript file <%s> to convert', ...
                    param_value)
            end
            gsData.psFile = param_value;
            
            % path to pdf file to create
        case 'pdffile'
            %verify we can create file at that location
            pdf_fid = fopen(param_value,'w');
            if pdf_fid < 0
                error('ps2pdf:invalidPDFFIle', ...
                    'Can not open <%s> for writing', ...
                    param_value);
            end
            fclose(pdf_fid);
            %delete temp file we created
            delete(param_value);
            gsData.pdfFile = param_value;
            
            % full path to gs executable
        case 'gscommand'
            if ~exist(param_value, 'file')
                error('ps2pdf:ghostscriptCommand', ...
                    'Can not find Ghostscript executable (''gscommand'') <%s>',...
                    param_value)
            end
            if ispc && ~isempty(findstr(param_value, ' '))
                param_value = ['"' param_value '"']; %#ok<AGROW>
            end
            gsData.cmd = param_value;
            
            % full path to gs font dir
        case 'gsfontpath'
            if ~exist(param_value, 'dir')
                error('ps2pdf:ghostscriptFontPath', ...
                    'Can not find the directory <%s> for Ghostscript fonts (''gsfontpath'')', ...
                    param_value)
            end
            gsData.fontPath = param_value;
            
            % full path to gs lib dir
        case 'gslibpath'
            if ~exist(param_value, 'dir')
                error('ps2pdf:ghostscriptLibPath', ...
                    'Can not find the directory <%s> for Ghostscript library files (''gslibpath'')', ...
                    param_value)
            end
            gsData.libPath = param_value;
            
            % paper size
        case 'gspapersize'
            idx = strcmpi(param_value, gsData.paperSizes);
            if ~any(idx)
                warning('ps2pdf:papersize', ...
                    '''gspapersize'' value <%s> not found in the list of known sizes, ignoring it.', param_value);
            else
                gsData.paperSize = gsData.paperSizes{idx};
            end
            
            % deletePSFile
        case 'deletepsfile'
            if isnumeric(param_value)
                gsData.deletePSFile = param_value;
            else
                warning('ps2pdf:deletepsfile', ...
                    '''deletepsfile'' value <%s> class <%s> should be numeric, defaulting to 0', ...
                    param_value, class(param_value));
            end
            
            % verbose
        case 'verbose'
            if isnumeric(param_value)
                gsData.verbose = param_value;
            else
                warning('ps2pdf:verbose', ...
                    '''verbose'' value <%s> class <%s> should be numeric, defaulting to 0', ...
                    param_value, class(param_value));
            end
            
        otherwise
            if isnumeric(param_value)
                param_value = num2str(param_value);
            end
            warning('ps2pdf:unknown', ...
                'ignoring unknown parameter <%s> with value <%s>.', param_arg, param_value);
    end
end

if ~isfield(gsData, 'psFile')
    error('ps2pdf:noInputFile', ...
        'No input (psfile) file specified');
end

if ~isfield(gsData, 'pdfFile')
    error('ps2pdf:noOutputFile', ...
        'No output (pdffile) file specified');
end

if ~isfield(gsData, 'cmd')
    if isdeployed
        error('ps2pdf:deployedNeedsGhostscript', ...
            'In order to use ''ps2pdf'' in a deployed application you must provide the path to a separate instance of Ghostscript.');
    end
    
    % updated code to find ghostscript - look for gs8x first,
    % then try old location. Don't depend on MATLAB version #
    ghostDir = fullfile( matlabroot, 'sys', 'gs8x' );
    if ~exist(ghostDir, 'dir')
        [gsCmd, ghostDir] = Local_GetOldGhostscript();
        gsData.cmd = gsCmd;
    else
        gsData.cmd = Local_GetGscriptFcnHandle;
        if ~isempty(gsData.cmd)
            gsData.useBuiltin = 1; % use builtin Ghostscript
        end
    end
    if ~exist(ghostDir, 'dir')
        error('ps2pdf:ghostscriptCommand', ...
            'Can not find Ghostscript installed with MATLAB in <%s>',...
            ghostDir);
    end
    
    if ~isempty(gsData.cmd)
        % if using MATLAB's version of GhostScript, use same set of fonts and library files
        if isfield(gsData, 'fontPath') || isfield(gsData, 'libPath')
            warning('ps2pdf:ghostscriptPathOverride', ...
                'Using MATLAB''s version of Ghostscript; overriding ''gsfontpath'' and ''gslibpath'' to use builtin MATLAB version');
        end
        gsData.fontPath = fullfile( ghostDir, 'fonts', '');
        gsData.libPath = fullfile( ghostDir, 'ps_files', '');
    else
        error('ps2pdf:noGhostscriptCommand', ...
            'Can not find Ghostscript program in MATLAB');
    end
else
    % if gscommandpath was specified,
    if ~isfield(gsData, 'fontPath') || ~isfield(gsData, 'libPath')
        warning('ps2pdf:ghostscriptCommandSuggestion', ...
            ['When specifying a Ghostscript executable (''gscommand'') you should also '...
            'specify both the ''gsfontpath'' and ''gslibpath'' locations']);
        
    end
end
end

%local function to create the input file needed for Ghostscript
function gsData = LocalCreateResponseFile(gsData)
% open a response file to write out Ghostscript commands
rsp_file = [tempname '.rsp'];
rsp_fid = fopen (rsp_file, 'w');

if (rsp_fid < 0)
    error('ps2pdf:responseFileCreate', 'Unable to create response file')
end

fprintf(rsp_fid, '-dBATCH -dNOPAUSE\n');
if ~gsData.verbose
    fprintf(rsp_fid, '-q\n');
end
if isfield(gsData, 'libPath')
    fprintf(rsp_fid, '-I"%s"\n', gsData.libPath);
end
if isfield(gsData, 'fontPath')
    fprintf(rsp_fid, '-I"%s"\n', gsData.fontPath);
end
if isfield(gsData, 'paperSize')
    fprintf( rsp_fid, '-sPAPERSIZE=%s\n', gsData.paperSize );
end
fprintf(rsp_fid, '-sOutputFile="%s"\n', gsData.pdfFile);
fprintf(rsp_fid, '-sDEVICE=%s\n', 'pdfwrite');
fclose(rsp_fid);
gsData.responseFile = rsp_file;
end

%local function to get a handle to MATLAB's Ghostscript implementation
%NOTE: this may change or be removed in future releases
function gs = Local_GetGscriptFcnHandle()
gs = '';
p = which('-all', 'gscript');
if ~isempty(p)
    p = p{1};
    fpath = fileparts(p);
    olddir = cd(fpath);
    gs = @gscript;
    cd(olddir);
end
end

% local function to try and get location of Ghostscript in older MATLAB
function [gsCmd, ghostDir] = Local_GetOldGhostscript
ghostDir = fullfile( matlabroot, 'sys', 'ghostscript' );
gsCmd = '';
if ispc
    if exist(fullfile(ghostDir,'bin','win32','gs.exe'), 'file')
        gsCmd = fullfile(ghostDir,'bin','win32','gs.exe');
    end
else
    if exist(fullfile(ghostDir,'bin',lower(computer),'gs'), 'file')
        gsCmd = fullfile(ghostDir,'bin',lower(computer), 'gs');
    end
end
gsCmd = ['"' gsCmd '"'];
end

function c = linewrap(s, maxchars)
%LINEWRAP Separate a single string into multiple strings
%   C = LINEWRAP(S, MAXCHARS) separates a single string into multiple
%   strings by separating the input string, S, on word breaks.  S must be a
%   single-row char array. MAXCHARS is a nonnegative integer scalar
%   specifying the maximum length of the broken string.  C is a cell array
%   of strings.
%
%   C = LINEWRAP(S) is the same as C = LINEWRAP(S, 80).
%
%   Note: Words longer than MAXCHARS are not broken into separate lines.
%   This means that C may contain strings longer than MAXCHARS.
%
%   This implementation was inspired a blog posting about a Java line
%   wrapping function:
%   http://joust.kano.net/weblog/archives/000060.html
%   In particular, the regular expression used here is the one mentioned in
%   Jeremy Stein's comment.
%
%   Example
%       s = 'Image courtesy of Joe and Frank Hardy, MIT, 1993.'
%       c = linewrap(s, 40)
%
%   See also TEXTWRAP.

% Steven L. Eddins
% $Revision: 1.7 $  $Date: 2006/02/08 16:54:51 $

error(nargchk(1, 2, nargin));

bad_s = ~ischar(s) || (ndims(s) > 2) || (size(s, 1) ~= 1);
if bad_s
    error('S must be a single-row char array.');
end

if nargin < 2
    % Default value for second input argument.
    maxchars = 80;
end

% Trim leading and trailing whitespace.
s = strtrim(s);

% Form the desired regular expression from maxchars.
exp = sprintf('(\\S\\S{%d,}|.{1,%d})(?:\\s+|$)', maxchars, maxchars);

% Interpretation of regular expression (for maxchars = 80):
% '(\\S\\S{80,}|.{1,80})(?:\\s+|$)'
%
% Match either a non-whitespace character followed by 80 or more
% non-whitespace characters, OR any sequence of between 1 and 80
% characters; all followed by either one or more whitespace characters OR
% end-of-line.

tokens = regexp(s, exp, 'tokens').';

% Each element if the cell array tokens is single-element cell array
% containing a string.  Convert this to a cell array of strings.
get_contents = @(f) f{1};
c = cellfun(get_contents, tokens, 'UniformOutput', false);

% Remove trailing whitespace characters from strings in c.  This can happen
% if multiple whitespace characters separated the last word on a line from
% the first word on the following line.
c = deblank(c);


end
