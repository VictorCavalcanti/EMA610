function [MAC] = mac(PHI1,PHI2)
MAC = ((PHI1'*PHI2).^2)./(diag(PHI1'*PHI1)*diag(PHI2'*PHI2)');
figure;
imagesc(MAC);
colormap(parula(32)); colorbar;
axis equal;
xlim([0.5, numel(PHI2(1,:))+0.5]);
ylim([0.5, numel(PHI1(:,1))+0.5]);
set(gca,'XTick',1:numel(PHI2(1,:)));
set(gca,'YTick',1:numel(PHI1(:,1)));
% ylabel('FEM Mode #');
% xlabel('Test Mode #');
% title('Modal Assurance Criterion (MAC) between Test and FEM modes');

end
