%code to get T location across multiple sets

sets = 400:514;
all_locations_by_cnd = cell(1,length(sets));

CNDFile = 'Z:\eblab\Cortex Programs\Contextual\TL2.cnd';
cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

cnd2itm = [];
for line = 1:size(cndfil)-1;
    str = textscan(cndfil(line+1,:),'%d');
    cnd2itm(line) = str{1}(2);
end
%%
for it = 1:length(sets)
    if exist(['Z:\eblab\Tloc\Tloc' num2str(sets(it)) '.xls'],'file')
        tloc = xlsread(['R:\Buffalo Lab\eblab\Tloc\Tloc' num2str(sets(it)) '.xls']); %arranged by item #
        
        locs = NaN(2,length(cnd2itm)); %[x;y]
        for cnd = 1:length(cnd2itm)
            locs(1,cnd)=tloc(cnd2itm(cnd),3);
            locs(2,cnd)=tloc(cnd2itm(cnd),4);
        end
        all_locations_by_cnd{it} = locs;
    else
        all_locations_by_cnd{it} = NaN;
    end
end

save('Tlocations','all_locations_by_cnd','sets')
    
