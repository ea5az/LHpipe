
%subplot(2,2,4)
%scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, (tab.rates > 0.6), 'filled')
% 
% obj = fitgmdist(X,3,'Replicates',20,'Options',options);
% ppb = posterior(obj , [tab.amps + 1 , tab.jitter , tab.rates]);
% [~,idx] = max(ppb');
% subplot(4,2,5)
% colormap(jet); scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, idx, 'filled')
% ylabel(sprintf('AIC: %4.0f , BIC: %4.0f',obj.AIC , obj.BIC))
% subplot(4,2,6)
% scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, (tab.rates > 0.1) + (tab.rates > 0.6), 'filled')
% 
% obj = fitgmdist(X,4,'Replicates',20,'Options',options);
% ppb = posterior(obj , [tab.amps + 1 , tab.jitter , tab.rates]);
% [~,idx] = max(ppb');
% subplot(4,2,7)
% colormap(jet); scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, idx, 'filled')
% ylabel(sprintf('AIC: %4.0f , BIC: %4.0f',obj.AIC , obj.BIC))
% subplot(4,2,8)
% scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, (tab.rates > 0.1) + (tab.rates > 0.6), 'filled')