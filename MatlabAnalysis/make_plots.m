function make_plots()
% Change default fonts.
set(0,'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesFontWeight', 'Light');
set(0,'DefaultAxesFontSize', 10);
set(0,'DefaultTextFontName', 'Arial');
set(0,'DefaultTextFontWeight', 'Light');
set(0,'DefaultTextFontSize', 10);
set(0,'DefaultAxesTickDir', 'out');
set(0,'DefaultAxesTickLength', [.01 0.1]);

set(0,'DefaultAxesBox', 'off');
set(0,'DefaultAxesXMinorTick', 'on');
set(0,'DefaultAxesYMinorTick', 'on');


LEGENDFONTSIZE=8;

% load stuff
results

kinds = {'indiv', 'multi'};
short_families = {'afro', 'austro', 'indo', 'niger', 'nilo', 'sino'};
long_families = {'Afro-Asiatic', 'Austronesian', 'Indo-European', 'Niger-Congo', 'Nilo-Saharan', 'Sino-Tibetan'};
methods = {'geographic', 'genetic', 'feature', 'combination'};
%styles = {'-', '--', '.', '-.', ':', '-:'};
%styles = {'-.ok', '-.sk', '-.*k', '-..k', '-.xk', '-.+k'};
styles = {'-r', '-b', '-g', '-c', '-m', '-y'};
for k=1:2
    kind = kinds{k};
    for i=1:4
        method = methods{i};
        
        % Do dynamics plots              
%         for j=1:6
%             %%% ADMIN
%             family = short_families{j};         
%             fig = figure();
%             
%             if k == 1
%                 key = strcat(method, '_', family);
%                 matrix = trans(key);
%             else
%                 key = strcat(method);
%                 matrix = trans(key);
%             end
%             
%             %%% HEATMAP
%             [x, y, width, height, aspectratio] = dynamics_plot_layout(1,1);
%             subplot('Position', [x, y, width, height]);
%             heatmap(matrix, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'},[],'ShowAllTicks',true);
%             set(gca,'ticklength',[0,0]);
%             axis square;
%             colormap(hot);
%             xlabel('Target word order');
%             ylabel('Source word order');
%             title('Word order transition probabilities');
%                         
%             %%% STABILITIES
%             [x, y, width, height, aspectratio] = dynamics_plot_layout(1,2);
%             subplot('Position', [x, y, width, height]);           
%             stab_samples = all_stabs(key);      
%             hold on        
%             for l=1:6
%                 stabs = stab_samples(:,l);
%                 % Invert things so higher values are more stable
%                 stabs = 1 ./ stabs;
%                 % Scale to years
%                 stabs = 10000 .* stabs;
%                 [y, x] = ksdensity(stabs);                
%                 plot(x,y, styles{l},'Linewidth',1);
%             end
%             hold off;
%             title('Posterior distribution of stability parameters');
%             ylabel('Probability');
%             xlabel('Mean time between changes (thousands of years)');
%             xlim([0 30000]);
%             ylim([0 26/10000]);
%             set(gca,'XTick',0:5000:30000)
%             set(gca,'XTickLabel',0:5:30)           
%             h_legend = legend('SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV');                    
%             set(h_legend,'FontSize',LEGENDFONTSIZE);
%         
%             %%% SAVE
%             filename = strcat('plots/', key, '_dynamics.png');
%             set(gcf,'PaperPositionMode','auto')
%             set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 10*aspectratio]);
%             print('-dpng',filename);
%             close(fig);
%          end
%         
        % Do ancestral plots
        
%         for prior=1:5
        prior = 4
            fig = figure();
            for j=1:6
                family = short_families{j};    

                key = strcat(method, '_', family);
                row = floor(j/4) + 1;
                col = j - floor(j/4)*3;
                [x, y, width, height, aspectratio] = ancestral_plot_layout(row,col);
                subplot('Position', [x, y, width, height]);
                switch prior
                    case 1
                        priortype = 'uniform';
                        if k == 1
                            vector = indiv_uniform_ancestrals(key);
                        else
                            vector = multi_uniform_ancestrals(key);
                        end
                    case 2
                        priortype = 'fuzzy';
                        if k == 1
                            vector = indiv_fuzzy_ancestrals(key);
                        else
                            vector = multi_fuzzy_ancestrals(key);
                        end
                    case 3
                        priortype = 'stationary';
                        if k == 1
                            vector = indiv_stationary_ancestrals(key);
                        else
                            vector = multi_stationary_ancestrals(key);
                        end
                    case 4
                        priortype = 'present';
                        if k == 1
                            vector = indiv_uniform_ancestrals(key);
                        else
                            vector = multi_uniform_ancestrals(key);
                        end
                        presentdist = [0.48, 0.41, 0.08, 0.02, 0.01, 0.0];
                        norm = 0;
                        for unusedindex=1:6
                            vector(unusedindex) = vector(unusedindex) * presentdist(unusedindex);
                            norm = norm + vector(unusedindex);
                        end
                        vector = vector ./ norm;
                    case 5
                        priortype = 'tom';
                        if k == 1
                            uniform = indiv_uni_ancestrals(key);
                        else
                            uniform = multi_uni_ancestrals(key);
                        end
                        vector = [0,0,0,0,0,0];
                        norm = 0;
                        for ui=1:6
                            denom = 0;
                            for foo=1:6
                                if foo ~= ui
                                    denom = denom + uniform(foo)*(1.0/5.0);
                                end
                            end
                            vector(ui) = uniform(ui) / denom;
                            %norm = norm + vector(ui);
                        end
                        vector = log(vector);
                        %vector = vector ./ norm;                                                
                end
                
                bar(vector, 'black');            
                box off
                if prior == 5
                    axis([0,7,-5,2]);
                    ylabel('Evidence (nats)');
                else
                    axis([0,7,0.0,1.0]);
                    ylabel('Posterior probability');
                end
                axis square;

                %set(gca,'xtick',[])
                %set(gca,'ytick',[])
                set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'});
                set(gca,'XMinorTick', 'off');
                set(gca,'YMinorTick', 'off');                
                %set(gca,'LooseInset',get(gca,'TightInset'))
                %set(gca,'Position',[.05 .15, .9 .7])
                title(long_families{j});      
            end
            filename = strcat('plots/', kind , '_', method, '_', priortype, '_ancestral.png');
            set(gcf,'PaperPositionMode','auto')
            %set(fig, 'Position', [0 0 100 15])
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 10*aspectratio]);
            print('-dpng',filename);
            close(fig);
%         end       

        
        
     end


%        Do long term plot
set(0,'DefaultAxesFontSize', 13);
fig = figure();
[x, y, width, height, aspectratio] = longterm_plot_layout(1,1);
subplot('Position', [x, y, width, height]);        
title('Posterior probability of common ancestor word orders');
xlabel('Age of common ancestor (thousands of years before present)');
ylabel('Probability');
axis([35,100,0.1,0.2]);    
hold on;
for i=1:6
    ages = common_ancestor_ages ./ 1000;
    data = common_ancestor_dists(:,i);            
    plot(ages(1:5:end), data(1:5:end), styles{i}, 'Linewidth',1);
end
hold off;
h_legend = legend('SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV','Location','SouthEast');
set(h_legend,'FontSize',LEGENDFONTSIZE+1);

[x, y, width, height, aspectratio] = longterm_plot_layout(2,1);
subplot('Position', [x, y, width, height]);        
title('Predicted evolution of present word order distribution');
xlabel('Time into future (thousands of years after present)');
ylabel('Probability');

axis([0,300,0,0.8]);
future = [0.48, 0.41, 0.08, 0.02, 0.01, 0];

q = Q('combination');
p = expm(100.0/10000.0 * q);
for i=1:300
    newdist = [0,0,0,0,0,0];
    for final=1:6
        for initial=1:6
            newdist(final) = newdist(final) + future(i,initial)*p(initial,final);
        end
    end
    future = [future ; newdist];
end
hold on;
for i=1:6
    % Christ, I can't stand this ridiculous shit!
    % I really have to assign the result of every slice operation
    % to an intermediate variable so that I can do more slicing?!
    ages = 0:300;
    data = future(:,i);            
    plot(ages(1:5:end), data(1:5:end), styles{i},'Linewidth',1);
end
hold off;
%legend('SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV');
filename = strcat('plots/common_ancestral.png');        
set(gcf,'PaperPositionMode','auto')
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 10*aspectratio]);
print('-dpng',filename);
close(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(0,'DefaultAxesFontSize', 17);
    for k=1:2
        kind = kinds{k};
        % Iterate over methods
        for m=1:4
            method = methods{m};
            fig = figure();
            cell = 1;
            % Iterate over families
            for j=1:6
                family = short_families{j};
                key = strcat(method, '_',family);

                row = j;
                col = 1;
                [x, y, width, height, aspectratio] = sliding_plot_layout(row,col);
                subplot('Position', [x, y, width, height]);
                if k==1
                    vector = indiv_uniform_ancestrals(key);
                else
                    vector = multi_uniform_ancestrals(key);
                end
                bar(compute_evidence(vector), 'black');
                box off;
                %xlim([0 7]);
                if k==1
                    axis([0,7,-1,1]);
                else
                    if m == 1
                        axis([0,7,-9,3]);
                    elseif m == 3
                        axis([0,7,-13,4]);
                    else
                        axis([0,7,-8,3]);
                    end
                end
                ylabel('Evidence (bits)','FontSize',17)%,'FontWeight','bold');
                %axis square;
                set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'},'FontSize',12);
                %set(gca,'xtick',[])
                %set(gca,'ytick',[])
                set(gca,'XMinorTick', 'off');
                set(gca,'YMinorTick', 'off');

                title(long_families{j},'FontSize',17)%,'FontWeight','bold');
                %h = gca;
                %rotateticklabel(gca,0);
                cell = cell + 1;

                row = j;
                col = 2;
                [x, y, width, height, aspectratio] = sliding_plot_layout(row,col);
                subplot('Position', [x, y, width, height]);
                hold on;
                if k==1
                    plot(indiv_sliding_priors(key),'Linewidth',1);
                else
                    plot(multi_sliding_priors(key),'Linewidth',1);
                end
                if j==1
                    h_legend = legend('SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV','Location','NorthEast');
                    set(h_legend,'FontSize',12);
                end
                %title(long_families{j});
                ylabel('Posterior probability')%,'FontWeight','bold');
                xlabel('Prior')%,'FontWeight','bold');
                axis([0,100,0.0,1.0]);
                set(gca,'XTick',[0,100]);
                set(gca,'XTickLabel',{'Uniform','Stationary'})%,'FontWeight','bold')
                set(gca,'YTick',[0.0,0.2,0.4,0.6,0.8,1.0]); 
                foo = get(gca,'YTickLabel');
                set(gca,'YTickLabel',foo,'FontSize',12);               
                hold off;
                cell = cell + 1;
            end

            filename = strcat('plots/',kind,'_',method,'_sliding_priors.png');
            set(gcf,'PaperPositionMode','auto');
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 10*aspectratio]);
            print('-dpng',filename);
            close(fig);
        end
    end
end

end

function [evidence] = compute_evidence(posteriors)

evidence = [0,0,0,0,0,0];
norm = 0;
for i=1:6
    denom = 0;
    for j=1:6
        if j ~= i
            denom = denom + posteriors(j)*(1.0/5.0);
        end
    end
    evidence(i) = posteriors(i) / denom;    
end
evidence = log2(evidence);

end
function [x, y, width, height, aspectratio] = dynamics_plot_layout(row, col)

unit2gap = 4.0;
aspectratio = (2+1*unit2gap) / (3+3*unit2gap);
horizspace = 1.0/(3+3*unit2gap);
vertspace = horizspace/aspectratio;
width = horizspace*unit2gap;
height = vertspace*unit2gap;
row=row-1;
if col == 1
    x = horizspace;
else
    x = 2*horizspace + width;
    width = 2*width;
end
y = vertspace + row*(height+vertspace);
    
end

function [x, y, width, height, aspectratio] = ancestral_plot_layout(row, col)

unit2gap = 4.0;
aspectratio = (3+2*unit2gap) / (4+3*unit2gap);
horizspace = 1.0/(4+3*unit2gap);
vertspace = horizspace/aspectratio;
width = horizspace*unit2gap;
height = vertspace*unit2gap;
x = horizspace + (col-1)*(horizspace+width);
y = vertspace + (row-1)*(vertspace+height);
    
end

function [x, y, width, height, aspectratio] = longterm_plot_layout(row, col)

unit2gap = 4.0;
aspectratio = (3+2*unit2gap) / (2+2*unit2gap);
horizspace = 1.0/(2+2*unit2gap);
vertspace = horizspace/aspectratio;
width = horizspace*unit2gap;
height = vertspace*unit2gap;
if row == 1
    y = 2*vertspace + height;
else
    y = vertspace;    
end
width = 2*width;
x = horizspace;
    
end

function [x, y, width, height, aspectratio] = sliding_plot_layout(row, col)

unit2gap = 3.75;
aspectratio = (7+6*unit2gap) / (3+3*unit2gap);
horizspace = 1.0/(3+3*unit2gap);
vertspace = horizspace/aspectratio;
width = horizspace*unit2gap;
height = vertspace*unit2gap;
row = 6-row;
if col == 1
    x = horizspace;
else
    x = 2*horizspace + width;
    width = 2*width;
end
y = vertspace + row*(height+vertspace);
    
end

