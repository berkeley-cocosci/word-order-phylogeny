% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Verdana')
set(0,'DefaultAxesFontSize', 10)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Verdana')
set(0,'DefaultTextFontSize', 10)

% load stuff
results

kinds = {'indiv', 'multi'};
short_families = {'afro', 'austro', 'indo', 'niger', 'nilo', 'sino'};
long_families = {'Afro-Asiatic', 'Austronesian', 'Indo-European', 'Niger-Congo', 'Nilo-Saharan', 'Sino-Tibetan'};
methods = {'geographic', 'genetic', 'feature', 'combination'};

for k=1:2
    kind = kinds{k};
    for i=1:4
        method = methods{i};
        
        % Do dynamics plots      
%         for j=1:6
%             family = short_families{j};
%           
%             fig = figure('units','normalized','outerposition',[0 0 1 1]);
%             if k == 1
%                 key = strcat(method, '_', family);
%                 matrix = trans(key);
%             else
%                 key = strcat(method);
%                 matrix = trans(key);
%             end
%             %subaxis(1,2,1,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
%             subaxis(1, 2, 1)
%             heatmap(matrix, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'});
%             colormap(flipud(gray));
%             axis square;
%             xlabel('Target word order');
%             ylabel('Source word order');
%             title('Word order transition probabilities');
%             stab_samples = all_stabs(key);       
%             %subaxis(1, 2, 2,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1)
%             subaxis(1, 2, 2)
%             hold on        
%             for l=1:6
%                 stabs = stab_samples(:,l)
%                 [y, x] = ksdensity(stabs);                
%                 plot(x,y);
%             end
%             hold off
%             filename = strcat('plots/', key, '_dynamics.png');        
%             print('-dpng',filename);
%             close(fig);
%         end
%         
%         % Do ancestral plot
%         fig = figure('units','normalized','outerposition',[0 0 1 1]);
%         for j=1:6
%             family = short_families{j};    
%                       
%             key = strcat(method, '_', family);
%             %row = floor(j/3) + 1
%             %col = j - floor(j/3)
%             subaxis(2, 3, j);
%             if k == 1
%                 vector = indiv_fuzzy_ancestrals(key);
%             else
%                 vector = multi_fuzzy_ancestrals(key);
%             end
%             bar(log(vector ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black');            
%             ylabel('Evidence');
%             axis([0,7,-1.2,0.8]);
%             axis square;
%             set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'});
%             %set(gca,'LooseInset',get(gca,'TightInset'))
%             %set(gca,'Position',[.05 .15, .9 .7])
%             title(long_families{j});       
%         end
%         filename = strcat('plots/', kind , '_', method, '_fuzzy_ancestral.png');        
%         print('-dpng',filename);
%         close(fig);
        
        
        % Do long term plot
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        subaxis(2, 1, 1);
        hold on;
        for i=1:6
            plot(common_ancestor_ages, common_ancestor_dists(:,i));
        end
        hold off;
        
        subaxis(2, 1, 2);
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
            plot(future(:,i));
        end
        hold off;

        filename = strcat('plots/common_ancestral.png');        
        print('-dpng',filename);
        close(fig);
        
        %         for j=1:6
%             family = short_families{j};         
%             fig = figure;
%             key = strcat(method, '_', family);
%             if k == 1                
%                 vector = indiv_fuzzy_stationary(key);
%                 filename = strcat('plots/indiv_', key, '_longterm.png');        
%             else                
%                 vector = multi_fuzzy_stationary(key);
%                 filename = strcat('plots/multi_', key, '_longterm.png');        
%             end
%             subaxis(1,2,1,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
%             
%             
%             title('Word order transition probabilities');           
%             print('-dpng',filename);
%             close(fig);
%         end        
        
        
    end
end
        
%heatmap(trans('combination_indo'), {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
%%%%%%%%%%%%%%%%%%%%%%%%5

% figure
% 
% heatmap(common_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% %print('-dpng','common_feature_transition.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%5
% 
% subaxis(2,4,2,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_afro_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% %set(gca,'Position',[.05 .15, .9 .7])
% title('Afro-Asiatic')
% %print('-dpng','common_afro_feature_postancestral.png')
% 
% subaxis(2,4,3,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_austro_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Austronesian')
% %print('-dpng','common_austro_feature_postancestral.png')
% 
% subaxis(2,4,4,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_indo_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Indo-European')
% %print('-dpng','common_indo_feature_postancestral.png')
% 
% subaxis(2,4,6,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_nilo_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Nilo-Saharan')
% %print('-dpng','common_nilo_feature_postancestral.png')
% 
% subaxis(2,4,7,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_niger_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Niger-Congo')
% %print('-dpng','common_niger_feature_postancestral.png')
% 
% subaxis(2,4,8,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_sino_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Sino-Tibetan')
% %print('-dpng','common_sino_feature_postancestral.png')
% 
% %%%%%%%%%%%%%%%%%
% 
% subaxis(2,4,5,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(common_feature_stationary, 'black')
% xlabel('Word order')
% ylabel('Stationary probability')
% axis([0,7,0.0,1.0])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Stationary word order probabilities')
% %print('-dpng','common_feature_stationary.png')
% print('-dpng','everything_feature.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% subaxis(2,4,1,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(common_family_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% %print('-dpng','common_family_transition.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%5
% 
% subaxis(2,4,2,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_afro_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% %set(gca,'Position',[.05 .15, .9 .7])
% title('Afro-Asiatic')
% %print('-dpng','common_afro_family_postancestral.png')
% 
% subaxis(2,4,3,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_austro_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Austronesian')
% %print('-dpng','common_austro_family_postancestral.png')
% 
% subaxis(2,4,4,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_indo_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Indo-European')
% %print('-dpng','common_indo_family_postancestral.png')
% 
% subaxis(2,4,6,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_nilo_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Nilo-Saharan')
% %print('-dpng','common_nilo_family_postancestral.png')
% 
% subaxis(2,4,7,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_niger_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Niger-Congo')
% %print('-dpng','common_niger_family_postancestral.png')
% 
% subaxis(2,4,8,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_sino_family_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Sino-Tibetan')
% %print('-dpng','common_sino_family_postancestral.png')
% 
% %%%%%%%%%%%%%%%%%
% 
% subaxis(2,4,5,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(common_family_stationary, 'black')
% xlabel('Word order')
% ylabel('Stationary probability')
% axis([0,7,0.0,1.0])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Stationary word order probabilities')
% %print('-dpng','common_family_stationary.png')
% print('-dpng','everything_family.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% subaxis(2,4,1,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(common_distance_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% %print('-dpng','common_distance_transition.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%5
% 
% subaxis(2,4,2,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_afro_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% %set(gca,'Position',[.05 .15, .9 .7])
% title('Afro-Asiatic')
% %print('-dpng','common_afro_distance_postancestral.png')
% 
% subaxis(2,4,3,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_austro_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Austronesian')
% %print('-dpng','common_austro_distance_postancestral.png')
% 
% subaxis(2,4,4,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_indo_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Indo-European')
% %print('-dpng','common_indo_distance_postancestral.png')
% 
% subaxis(2,4,6,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_nilo_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Nilo-Saharan')
% %print('-dpng','common_nilo_distance_postancestral.png')
% 
% subaxis(2,4,7,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_niger_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% %axis tight
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Niger-Congo')
% %print('-dpng','common_niger_distance_postancestral.png')
% 
% subaxis(2,4,8,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(common_sino_distance_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% %xlabel('Word order')
% %ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% %set(gca,'LooseInset',get(gca,'TightInset'))
% title('Sino-Tibetan')
% %print('-dpng','common_sino_distance_postancestral.png')
% 
% %%%%%%%%%%%%%%%%%
% 
% subaxis(2,4,5,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(common_distance_stationary, 'black')
% xlabel('Word order')
% ylabel('Stationary probability')
% axis([0,7,0.0,1.0])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Stationary word order probabilities')
% %print('-dpng','common_distance_stationary.png')
% print('-dpng','everything_distance.png')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% subaxis(3,4,1,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(afro_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%5
% 
% subaxis(3,4,2,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(afro_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Afro-Asiatic')
% 
% subaxis(3,4,3,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(austro_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% subaxis(3,4,4,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(austro_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Austronesian')
% 
% subaxis(3,4,5,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(indo_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% subaxis(3,4,6,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(indo_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Indo-European')
% 
% subaxis(3,4,7,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(nilo_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% subaxis(3,4,8,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(nilo_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Nilo-Saharan')
% 
% subaxis(3,4,9,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(niger_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% subaxis(3,4,10,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(niger_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Niger-Congo')
% 
% subaxis(3,4,11,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% heatmap(sino_feature_transition, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'}, {'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% colormap(flipud(gray))
% axis square
% xlabel('Target word order')
% ylabel('Source word order')
% title('Word order transition probabilities')
% 
% subaxis(3,4,12,'Spacing', 0.03, 'Padding', 0, 'SpacingVert', 0.1);
% bar(log(sino_feature_ancestral ./ [1/6,1/6,1/6,1/6,1/6,1/6]), 'black')
% ylabel('Evidence')
% axis([0,7,-1.2,0.8])
% axis square
% set(gca,'XTickLabel',{'SOV', 'SVO', 'VSO', 'VOS', 'OVS', 'OSV'})
% title('Sino-Tibetan')
% print('-dpng','six_feature.png')