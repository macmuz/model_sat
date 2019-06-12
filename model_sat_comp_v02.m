clear all
close all

regions = {'Bothnian_Bay','Bothnian_Sea','Gulf_of_Finland','total'};
formatSpec = 'box%02d';

cat_folder = 'wyniki';
model_box_folder = 'wyniki/test_albedo';
sat_box_folder = 'wyniki/test_albedo';
model_reg_folder = 'wyniki/MOSAIC';
mosaic_folder = 'wyniki/MOSAIC';
thick_folder = 'wyniki/THICK-L4';
output_hi = 'wyniki/HI';
output_bins = 'wyniki/bins';

%cell with cathegories
cat_data=importdata([cat_folder '/cat.txt']);
cat = {[num2str(cat_data(2)), '-', num2str(cat_data(5)), ' cm']};

cat{end+1} = [num2str(cat_data(5)), '-', num2str(cat_data(6)), ' cm'];
cat{end+1} = [num2str(cat_data(6)), '-', num2str(cat_data(7)), ' cm'];
cat{end+1} = [num2str(cat_data(7)), '-', num2str(cat_data(8)), ' cm'];
% for i = 3:length(cat_data)-1
%     cat{end+1} = [num2str(cat_data(i)), '-', num2str(cat_data(i+1)), ' cm'];
% end
cat{end+1} = ['>', num2str(cat_data(end)), ' cm'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iteration through boxes
for i = 1:13
%     break;
    nbox = sprintf(formatSpec,i);
    
    model_box_data = importdata([model_box_folder '/model_' nbox '.txt']);
    model_lvl_box_data = importdata([model_box_folder '/model_lvl_' nbox '.txt']);
    sat_box_data = importdata([sat_box_folder '/sat_' nbox '.txt']);
    mosaic_box_data = importdata([mosaic_folder '/sat_' nbox '.txt']);
    thick_box_data = importdata([thick_folder '/sat_' nbox '.txt']);
    
    model_box_bins_data = importdata([model_box_folder '/model_bins_' nbox '.txt']);
    model_lvl_box_bins_data = importdata([model_box_folder '/model_lvl_bins_' nbox '.txt']);
    sat_box_bins_data = importdata([sat_box_folder '/sat_bins_' nbox '.txt']);
    mosaic_box_bins_data = importdata([mosaic_folder '/sat_bins_' nbox '.txt']);
    thick_box_bins_data = importdata([thick_folder '/sat_bins_' nbox '.txt']);
    
    N1=max(size(model_box_data));
    model_time=datenum(model_box_data(1:N1,1),model_box_data(1:N1,2),model_box_data(1:N1,3));
    N2=max(size(sat_box_bins_data));
    sat_time=datenum(sat_box_data(1:N2,1),sat_box_data(1:N2,2),sat_box_data(1:N2,3),...
        sat_box_data(1:N2,4),sat_box_data(1:N2,5),sat_box_data(1:N2,6));
    N3=max(size(mosaic_box_data));
    mosaic_time=datenum(mosaic_box_data(1:N3,1),mosaic_box_data(1:N3,2),mosaic_box_data(1:N3,3),...
        mosaic_box_data(1:N3,4),mosaic_box_data(1:N3,5),mosaic_box_data(1:N3,6));
    N4=max(size(thick_box_data));
    thick_time=datenum(thick_box_data(1:N4,1),thick_box_data(1:N4,2),thick_box_data(1:N4,3),...
        thick_box_data(1:N4,4),thick_box_data(1:N4,5),thick_box_data(1:N4,6));

    model_hi = model_box_data(1:N1,4);
    model_lvl_hi = model_lvl_box_data(1:N1,4);
    sat_hi = sat_box_data(1:N2,7);
    mosaic_hi = mosaic_box_data(1:N3,7);
    thick_hi = thick_box_data(1:N4,7);
    
    model_bins=model_box_bins_data(1:N1,4:end)*100;
    model_lvl_bins=model_lvl_box_bins_data(1:N1,4:end)*100;
    sat_bins=sat_box_bins_data(1:N2,7:end)*100;
    mosaic_bins=mosaic_box_bins_data(1:N3,7:end)*100;
    thick_bins=thick_box_bins_data(1:N4,7:end)*100;
    
    model_bins = [sum(model_bins(1:N1,1:4),2) ,model_bins(1:N1,5:end)];
    model_lvl_bins = [sum(model_lvl_bins(1:N1,1:4),2),model_lvl_bins(1:N1,5:end)];
    mosaic_bins = [sum(mosaic_bins(1:N3,1:4),2) ,mosaic_bins(1:N3,5:end)];
    thick_bins = [sum(thick_bins(1:N4,1:4),2) ,thick_bins(1:N4,5:end)];
    sat_bins = [sum(sat_bins(1:N2,1:4),2),sat_bins(1:N2,5:end)];
    
    [~, ind] = unique(mosaic_time(:,1), 'rows');
    mosaic_hi=mosaic_hi(ind,:);
    mosaic_bins = mosaic_bins(ind,:);
    mosaic_time = mosaic_time(ind,:);
    [~, ind] = unique(thick_time(:,1), 'rows');
    thick_hi=thick_hi(ind,:);
    thick_bins = thick_bins(ind,:);
    thick_time = thick_time(ind,:);
    
    figure('Renderer', 'painters', 'Position', [10 300 900 600]);
    
    MyLineWidth=3;
    figure(1)
    set(gcf,'color','w');
    
    hold on;
    if ismember(i,[6,7,8,10,11])
        axis([datenum(2017,11,15) datenum(2018,6,15) 0 50]);
    else
        axis([datenum(2017,11,15) datenum(2018,6,15) 0 100]);
    end
    plot(model_time,model_lvl_hi,'-r','LineWidth',MyLineWidth);
    plot(sat_time,sat_hi,'+b','MarkerSize',MyLineWidth*3,'LineWidth',MyLineWidth);
    plot(mosaic_time,mosaic_hi,'+g','MarkerSize',MyLineWidth*3,'LineWidth',MyLineWidth);
    plot(thick_time,thick_hi,'+k','MarkerSize',MyLineWidth*3,'LineWidth',MyLineWidth);
    
    legend(['model lvl hi ',nbox],['sat ',nbox],...
        ['MOSAIC ',nbox],['THICK-L4 ',nbox],'Location','NorthEast');
    ylabel('mean ice thickness [cm]')
    set(gca,'FontSize',18);
    datetick('x','mmm');
    hold off;
    saveas(gcf, [ output_hi,'/hi_',nbox],'png')
    
    figure('Renderer', 'painters', 'Position', [920 300 900 600]);
    
    figure(2)
    colormap(jet);
    set(gcf,'color','w');

    barWidth=1;
    sp1=subplot(4,1,1);
    tmp=get(sp1,'Position');
    tmp(3) = 0.65;
    set(sp1,'Position',tmp);
    bar(model_time,model_lvl_bins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['model lvl ',nbox]);
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    set(gca,'FontSize',16,'XTick',[]);
    
    sp2=subplot(4,1,2);
    oldTime=sat_time;
    newTime=(sat_time(1):abs((sat_time(1)-sat_time(end))/180):sat_time(end));
    oldBins=sat_bins;
    for ii=1:min(size(oldBins)),
        cs=spline(oldTime',[oldBins(1,ii) oldBins(:,ii)' oldBins(end,ii)]);
        newBins(:,ii)=ppval(cs,newTime');
    end;
    tmp=get(sp2,'Position');
    tmp(3) = 0.65;
    set(sp2,'Position',tmp);
    bar(newTime,newBins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['sat ',nbox]);
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    set(gca,'FontSize',16,'XTick',[]);
    
    sp3=subplot(4,1,3);
    tmp=get(sp3,'Position');
    tmp(3) = 0.65;
    set(sp3,'Position',tmp);
    oldTime=mosaic_time;
    newTime=(mosaic_time(1):abs((mosaic_time(1)-mosaic_time(end))/180):mosaic_time(end));
    oldBins=mosaic_bins;
    for ii=1:min(size(oldBins)),
        cs=spline(oldTime',[oldBins(1,ii) oldBins(:,ii)' oldBins(end,ii)]);
        newBins(:,ii)=ppval(cs,newTime');
    end;
    
    bar(newTime,newBins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['mosaic ',nbox]);
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    set(gca,'FontSize',16,'XTick',[]);
    
    sp4=subplot(4,1,4);
    tmp=get(sp4,'Position');
    tmp(3) = 0.65;
    set(sp4,'Position',tmp);
    oldTime=thick_time;
    newTime=(thick_time(1):abs((thick_time(1)-thick_time(end))/180):thick_time(end));
    oldBins=thick_bins;
    for ii=1:min(size(oldBins)),
        cs=spline(oldTime',[oldBins(1,ii) oldBins(:,ii)' oldBins(end,ii)]);
        newBins(:,ii)=ppval(cs,newTime');
    end;
    
    bar(newTime,newBins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['thick-l4 ',nbox]);
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
%     yl=ylabel('percentage [%]','Position',[7.3700e+05 250.0000 -1]);
    yl=ylabel('percentage [%]','Units', 'Normalized','Position',[-0.1 2.6 0]);
    set(gca,'FontSize',16);
    lgd2=legend(cat,'Position',[0.8056 0.4400 0.1411 0.1600]);
    saveas(gcf, [ output_bins,'/bins_',nbox],'png')
    
 %return   
     close all

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iteration through regions
for i = 1:4
    nbox=regions{i};
    
    model_reg_data = importdata([model_reg_folder '/model_' nbox '.txt']);
    model_lvl_reg_data = importdata([model_reg_folder '/model_lvl_' nbox '.txt']);
    mosaic_reg_data = importdata([mosaic_folder '/sat_' nbox '.txt']);
    thick_reg_data = importdata([thick_folder '/sat_' nbox '.txt']);
    
    model_reg_bins_data = importdata([model_reg_folder '/model_bins_' nbox '.txt']);
    model_lvl_reg_bins_data = importdata([model_reg_folder '/model_lvl_bins_' nbox '.txt']);
    mosaic_reg_bins_data = importdata([mosaic_folder '/sat_bins_' nbox '.txt']);
    thick_reg_bins_data = importdata([thick_folder '/sat_bins_' nbox '.txt']);
    
    N1=max(size(model_reg_data));
    model_time=datenum(model_reg_data(1:N1,1),model_reg_data(1:N1,2),model_reg_data(1:N1,3));
    N3=max(size(mosaic_reg_data));
    mosaic_time=datenum(mosaic_reg_data(1:N3,1),mosaic_reg_data(1:N3,2),mosaic_reg_data(1:N3,3),...
        mosaic_reg_data(1:N3,4),mosaic_reg_data(1:N3,5),mosaic_reg_data(1:N3,6));
    N4=max(size(thick_reg_data));
    thick_time=datenum(thick_reg_data(1:N4,1),thick_reg_data(1:N4,2),thick_reg_data(1:N4,3),...
        thick_reg_data(1:N4,4),thick_reg_data(1:N4,5),thick_reg_data(1:N4,6));
    
    model_hi = model_reg_data(1:N1,4);
    model_lvl_hi = model_lvl_reg_data(1:N1,4);
    mosaic_hi = mosaic_reg_data(1:N3,7);
    thick_hi = thick_reg_data(1:N4,7);
    
    model_bins=model_reg_bins_data(1:N1,4:end)*100;
    model_lvl_bins=model_lvl_reg_bins_data(1:N1,4:end)*100;
    mosaic_bins=mosaic_reg_bins_data(1:N3,7:end)*100;
    thick_bins=thick_reg_bins_data(1:N4,7:end)*100;
    
    model_lvl_bins = [sum(model_lvl_bins(1:N1,1:4),2),model_lvl_bins(1:N1,5:end)];
    mosaic_bins = [sum(mosaic_bins(1:N3,1:4),2) ,mosaic_bins(1:N3,5:end)];
    thick_bins = [sum(thick_bins(1:N4,1:4),2) ,thick_bins(1:N4,5:end)];
    
    [~, ind] = unique(mosaic_time(:,1), 'rows');
    mosaic_hi=mosaic_hi(ind,:);
    mosaic_bins = mosaic_bins(ind,:);
    mosaic_time = mosaic_time(ind,:);
    [~, ind] = unique(thick_time(:,1), 'rows');
    thick_hi=thick_hi(ind,:);
    thick_bins = thick_bins(ind,:);
    thick_time = thick_time(ind,:);
    
    figure('Renderer', 'painters', 'Position', [10 300 900 600]);
    
    MyLineWidth=3;
    figure(1)
    set(gcf,'color','w');
    
    hold on;
    if i==1
        axis([datenum(2017,11,15) datenum(2018,6,15) 0 100]);
    else
        axis([datenum(2017,11,15) datenum(2018,6,15) 0 50]);
    end
    plot(model_time,model_lvl_hi,'-r','LineWidth',MyLineWidth);
    plot(mosaic_time,mosaic_hi,'+g','MarkerSize',MyLineWidth*3,'LineWidth',MyLineWidth);
    plot(thick_time,thick_hi,'+k','MarkerSize',MyLineWidth*3,'LineWidth',MyLineWidth);
    
    legend({['model lvl hi ',nbox],['MOSAIC ',nbox],...
        ['THICK-L4 ',nbox]},'Interpreter', 'none','Location','NorthEast');
    ylabel('mean ice thickness [cm]')
    set(gca,'FontSize',18);
    datetick('x','mmm');
    hold off;
    saveas(gcf, [ output_hi,'/hi_',nbox],'png')
    
    figure('Renderer', 'painters', 'Position', [920 300 900 600]);
    
    figure(2)
    colormap(jet);
    set(gcf,'color','w');

    barWidth=1;
    sp1=subplot(3,1,1);
    tmp=get(sp1,'Position');
    tmp(3) = 0.65;
    set(sp1,'Position',tmp);
    bar(model_time,model_lvl_bins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['model lvl ',nbox],'Interpreter', 'none');
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    set(gca,'FontSize',16,'XTick',[]);
    
    sp2=subplot(3,1,2);
    tmp=get(sp2,'Position');
    tmp(3) = 0.65;
    set(sp2,'Position',tmp);
    oldTime=mosaic_time;
    newTime=(mosaic_time(1):abs((mosaic_time(1)-mosaic_time(end))/180):mosaic_time(end));
    oldBins=mosaic_bins;
    for ii=1:min(size(oldBins)),
        cs=spline(oldTime',[oldBins(1,ii) oldBins(:,ii)' oldBins(end,ii)]);
        newBins(:,ii)=ppval(cs,newTime');
    end;
    
    bar(newTime,newBins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['mosaic ',nbox],'Interpreter', 'none');
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    set(gca,'FontSize',16,'XTick',[]);

    
    sp3=subplot(3,1,3);
    tmp=get(sp3,'Position');
    tmp(3) = 0.65;
    set(sp3,'Position',tmp);
    oldTime=thick_time;
    newTime=(thick_time(1):abs((thick_time(1)-thick_time(end))/180):thick_time(end));
    oldBins=thick_bins;
    for ii=1:min(size(oldBins)),
        cs=spline(oldTime',[oldBins(1,ii) oldBins(:,ii)' oldBins(end,ii)]);
        newBins(:,ii)=ppval(cs,newTime');
    end;
    
    bar(newTime,newBins,'stacked','BarWidth',barWidth);
    datetick('x','mmm');
    title(['thick-l4 ',nbox],'Interpreter', 'none');
    ylim([0 100]);
    xlim([datenum(2017,11,15) datenum(2018,6,15)]);
    
%     yl=ylabel('percentage [%]','Position',[7.3700e+05 250.0000 -1]);
    yl=ylabel('percentage [%]','Units', 'Normalized','Position',[-0.1 1.8 0]);
    set(gca,'FontSize',16);
    lgd2=legend(cat,'Position',[0.8056 0.4400 0.1411 0.1600]);
    saveas(gcf, [ output_bins,'/bins_',nbox],'png')
    
    
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

