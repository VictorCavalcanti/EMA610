function [CO]=corl8(PHI1,PHI2,M)
CO = PHI1'*M*PHI2;
figure;
imagesc(CO);
colormap(parula(32)); colorbar;
axis equal;
xlim([0.5, numel(PHI2(1,:))+0.5]);
ylim([0.5, numel(PHI1(:,1))+0.5]);
set(gca,'XTick',1:numel(PHI2(1,:)));
set(gca,'YTick',1:numel(PHI1(:,1)));

end