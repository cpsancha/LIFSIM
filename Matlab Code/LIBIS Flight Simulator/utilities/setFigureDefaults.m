function []= setFigureDefaults()
%Llamar a la función al final de la figura para tener el mismo formato en
%todas. Copiar el codigo comentado abajo para guardar la figura
%automaticamente en un archivo al que llamaremos desde al Lyx. Así, si se
%cambia algo en el matlab, el lyx lo substituira automaticamente por la
%nueva version.
     grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     grafWidth   = 16;
     grafAR      = 0.7;
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
%      saveas(gcf,format_Grafico,'epsc'); 
end

%% Codigo para guardar automaticamente. Copiar para cada figura

%  format_Grafico = 'nombre con el que quiero guardar la figura';
%  format_Grafico = [pwd filesep 'Carpeta donde quiero guardar la figura' filesep format_Grafico];
%  saveas(gcf,format_Grafico,'epsc'); 
