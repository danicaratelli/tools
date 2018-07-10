function IRFplotAllmod(Par,varargin)

Var_MP   = Par.Var_MP;      % index of the monetary policy variable
lag_plot = Par.lag_plot;    % maximal lag for plot

Var_plot = Par.Var_plotAll; % Variables of interest

level1 = Par.level1;        % confidence levels for the impulse responses
level2 = Par.level2;


shade1 = 0.85*ones(1,3);
shade2 = 0.7*ones(1,3);

fpatt = [0:lag_plot lag_plot:-1:0]';

for i = 1:length(varargin)
    irf = varargin{i};

    IRF       = irf.IRF;
    IRFG      = irf.IRFG;
    VarList   = irf.VarIdx;
    seriesVAR = irf.VarNames;

    VMP = zeros(max(VarList)); VMP(Var_MP)   = 1;   %% indicator for monetary policy variable
    VMP = VMP(VarList);
    VMP = find(VMP==1);

    VPl = zeros(max(VarList)); VPl(Var_plot) = 1;   %% indicator for variables of interest
    VPl = VPl(VarList);
    VPl = find(VPl==1);
    seriesPl = seriesVAR(VPl);

    IRF  = squeeze(IRF(VPl,VMP,1:lag_plot+1));
    IRFG = squeeze(IRFG(VPl,VMP,1:lag_plot+1,:));


    %==========================================================================
    % Confidence bands
    %==========================================================================

    m   = size(IRFG,3)-1;
    bb1 = ceil(m*(1-level1)/2);
    bb2 = ceil(m*(1-level2)/2);
    
    clear u_bound1 l_bound1 u_bound2 l_bound2
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

    for j = 1:length(seriesPl)
        
        axes1 = subplot(length(seriesPl),length(varargin),(j-1)*length(varargin)+i);
        set(axes1,'XTick',(0:12:lag_plot));
        set(axes1,'Xlim',[0 lag_plot])
%         set(axes1,'Ylim',Par.ylims(j,:));
        if i == 1 ylabel(axes1,seriesPl{j}); end
        box(axes1,'on');
        hold(axes1,'all')
        
        fpatt(:,2) = [u_bound1(j,:)'; flipud(l_bound1(j,:)')];
        fpatt(:,3) = [u_bound2(j,:)'; flipud(l_bound2(j,:)')];
        
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade1-0.1)
        patch(fpatt(:,1),fpatt(:,3),shade2,'EdgeColor',shade1-0.1)
        plot(0:lag_plot,IRF(j,:),'k.-',0:lag_plot,zeros(lag_plot+1),'k')
        if j==1 title(Par.Models{i});end
    end

end
% legend({['CI ',num2str(level1)],['CI ',num2str(level2)],'IRF'},'Orientation','Horizontal')
legend({[num2str(level1)],[num2str(level2)],'IRF'},'Orientation','Horizontal')
