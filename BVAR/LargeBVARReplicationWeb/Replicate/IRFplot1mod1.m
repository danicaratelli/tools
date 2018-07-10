function IRFplot1mod(Par,irf)

% Keys for the IRF plots (MEDIUM)
Var_plot = [33 115 113 2 3 6 20 25 51 109 125  129 87 72 73 77 78 83 93 104];

% configuration of the subplots
spDim = [5,4];

IRF       = irf.IRF; 
IRFG      = irf.IRFG; 
seriesVAR = irf.VarNames;
VarList   = irf.VarIdx;

Var_MP    = Par.Var_MP;     % index of the monetary policy variable
lag_plot  = Par.lag_plot;   % maximal lag for plot

level1 = Par.level1;        % confidence levels for the impulse responses
level2 = Par.level2;

VMP = zeros(max(Var_plot)); VMP(Var_MP)   = 1;   %% indicator for monetary policy variable
VMP = VMP(VarList);
VMP = find(VMP==1);
seriesMP = seriesVAR(VMP);

VPl = zeros(max(Var_plot)); VPl(Var_plot) = 1;   %% indicator for variables of interest
VPl = VPl(VarList);
VPl = find(VPl==1);
seriesPl = seriesVAR(VPl);


%==========================================================================
% Confidence bands
%==========================================================================
IRF  = squeeze(IRF(VPl,VMP,1:lag_plot+1));
IRFG = squeeze(IRFG(VPl,VMP,1:lag_plot+1,:));

m   = size(IRFG,3)-1;
bb1 = ceil(m*(1-level1)/2);
bb2 = ceil(m*(1-level2)/2);

for j = 1:size(IRFG,1)
    for l = 1:lag_plot+1
        temp = IRFG(j,l,:);
        temp = sort(temp);
        u_bound1(j,l) = temp(end-bb1);
        l_bound1(j,l) = temp(bb1);
        u_bound2(j,l) = temp(end-bb2);
        l_bound2(j,l) = temp(bb2);
    end;
end

%==========================================================================
% Plotting
%==========================================================================
seriesPl = changeMnem(seriesPl);

shade1 = 0.85*ones(1,3);
shade2 = 0.7*ones(1,3);

fpatt = [0:lag_plot lag_plot:-1:0]';
for i = 1:spDim(1)
    for j = 1:spDim(2)
        k = (i-1)*spDim(2)+j;
        
        axes1 = subplot(spDim(1),spDim(2),k); 
        set(axes1,'XTick',(0:12:lag_plot));
        set(axes1,'Xlim',[0 lag_plot])
        box(axes1,'on');
        hold(axes1,'all')

        fpatt(:,2) = [u_bound1(k,:)'; flipud(l_bound1(k,:)')];
        fpatt(:,3) = [u_bound2(k,:)'; flipud(l_bound2(k,:)')];
        
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade1-0.1)
        patch(fpatt(:,1),fpatt(:,3),shade2,'EdgeColor',shade1-0.1)
        plot(0:lag_plot,IRF(k,:),'k.-',0:lag_plot,zeros(lag_plot+1),'k')
        title(seriesPl{k},'FontSize',16)
    end
end
% legend({['CI ',num2str(level1)],['CI ',num2str(level2)],'IRF'},'Orientation','Horizontal')
legend({[num2str(level1)],[num2str(level2)],'IRF'},'Orientation','Horizontal')

