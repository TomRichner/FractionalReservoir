function out = fig_main8_psd(source_fig_root)
% FIG_MAIN8_PSD Restyle the selected archived PSD without re-estimating it.
% The selected mu7 archive has no native FIG or numerical PSD. Native FIGs
% from other runs have different artwork/limits and are not substitutes.
source=fullfile(source_fig_root,'fig_stim_engages_adaptation','bursting_psd.png');
assert(isfile(source),'fig_main8_psd:MissingSource','Missing selected PSD archive.');
fid=fopen(source,'rb'); guard=onCleanup(@()fclose(fid)); bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(bytes,'int8'));
digest=reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]); clear guard
assert(strcmpi(digest,'ad1408e6f780d41003c3e0c5c37f810bbc198ab48661b3f3f38a4674e8b72165'), ...
    'fig_main8_psd:UncalibratedRaster','Calibration requires the verified mu7 PSD archive.');
I=imread(source); cols=466:3463; rows=67:2137; C=I(rows,cols,:);
% Remove neutral old tick/text pixels. The legend lies wholly to the left of
% the first data frequency (0.3 Hz, source column approximately 940).
neutral=max(C,[],3)-min(C,[],3)<5;
for ch=1:3, plane=C(:,:,ch); plane(neutral)=255; C(:,:,ch)=plane; end
C(:,cols<900,:)=255;
fig=figure('Visible','off','Color','w','Position',[40 40 1300 760]);
ax=axes(fig,'Position',[.105 .16 .38 .74],'Tag','psd_model'); hold(ax,'on');
% A raster's spacing is uniform in log10 coordinates. Draw in those exact
% coordinates with native exponent ticks; no digitized/interpolated PSD.
image(ax,'XData',-1+(cols([1 end])-462)/3004*3, ...
    'YData',-(rows([1 end])-64)/2077.5*12,'CData',C,'Tag','archived_psd_artwork');
set(ax,'XLim',[-1 2],'YLim',[-12 0],'YDir','normal', ...
    'XTick',-1:2,'XTickLabel',{'10^{-1}','10^{0}','10^{1}','10^{2}'}, ...
    'YTick',[-10 -5 0],'YTickLabel',{'10^{-10}','10^{-5}','10^{0}'}, ...
    'FontSize',14,'LineWidth',1,'Box','off');
xlabel(ax,'frequency (Hz)','FontSize',14);
ylabel(ax,'Power spectral density of dendritic potential, x','FontSize',14,'Interpreter','none');
cmap=parula(6);
h1=plot(ax,NaN,NaN,'Color',cmap(1,:),'LineWidth',2);
h2=plot(ax,NaN,NaN,'Color',cmap(5,:),'LineWidth',2);
legend(ax,[h1 h2],{'no-stim','stim'},'Location','southwest','Box','off','FontSize',14);
bx=axes(fig,'Position',[.58 .16 .38 .74],'Tag','psd_human_empty', ...
    'XTick',[],'YTick',[],'XLim',[0 1],'YLim',[0 1], ...
    'Box','on','LineWidth',1,'FontSize',14);
for a=[ax bx]
    label='(A)'; if a==bx, label='(B)'; end
    text(a,-.10,1.06,label,'Units','normalized','FontSize',14, ...
        'Clipping','off','VerticalAlignment','bottom','Tag','panel_label');
end
notes={['Model PSD retains the selected mu7 archived trace pixels. No native FIG or numerical PSD ' ...
    'was saved for this source; other-run FIGs are not used. Data curves remain raster; axes, labels and legend are native.'], ...
    ['Raster calibration: 3551-by-2491 PNG; x pixels 462 to 3466 map to log10 frequency -1 to 2; ' ...
    'y pixels 64 to 2141.5 map to log10 PSD 0 to -12. Precision is about one source pixel. ' ...
    'Native axes use log10 coordinates and exponent tick labels to preserve uniform raster spacing exactly. ' ...
    'Neutral annotation pixels and the old legend left of the first data frequency are removed; scientific curve pixels are unchanged.'], ...
    'Panel B is an intentionally empty box. No patient data are shown. No titles; labels (A)/(B), 14-point fonts and axes linewidth 1.0.'};
out=struct('figs',fig,'source',{{source}},'notes',{notes});
end
