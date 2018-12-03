function [DATA, totDATA] = internal_norm_calista(DATA,data_type,totDATA)
DATA.original_data=totDATA;
DATA.min_totDATA=min(min(totDATA));
DATA.log_max_mRNA=log2(200);

switch data_type
    case 1
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        totDATA(totDATA>28)=28; %ct max
        DATA.ctmax=max(max(totDATA));
        log2Ex=DATA.ctmax-totDATA;
        DATA.max_totDATA=max(max(log2Ex));
        base=2^(DATA.log_max_mRNA/DATA.max_totDATA);
        totDATA=round(base.^log2Ex)-1;
    case 2
        totDATA=log2(totDATA+1);
        totDATA(totDATA>28)=28; %ct max
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        log2Ex=totDATA;
        DATA.max_totDATA=repmat(max(log2Ex),size(log2Ex,1),1);
        exponent=DATA.log_max_mRNA*(log2Ex./DATA.max_totDATA);
        totDATA=round(2.^exponent)-1;
    case 3
        totDATA=log2(totDATA+1);
        totDATA(totDATA>28)=28; %ct max
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        log2Ex=totDATA;
        DATA.max_totDATA=repmat(max(log2Ex),size(log2Ex,1),1);
        exponent=DATA.log_max_mRNA*(log2Ex./DATA.max_totDATA);
        totDATA=round(2.^exponent)-1;
    case 4
        totDATA=log2(log2(totDATA+1)+1);
        totDATA(totDATA>28)=28; %ct max
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        log2Ex=totDATA;
        DATA.max_totDATA=repmat(max(log2Ex),size(log2Ex,1),1);
        exponent=DATA.log_max_mRNA*(log2Ex./DATA.max_totDATA);
        totDATA=round(2.^exponent)-1;
    case 5
        totDATA=log2(log2(totDATA+1)+1);
        totDATA(totDATA>28)=28; %ct max
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        log2Ex=totDATA;
        DATA.max_totDATA=repmat(max(log2Ex),size(log2Ex,1),1);
        exponent=DATA.log_max_mRNA*(log2Ex./DATA.max_totDATA);
        totDATA=round(2.^exponent)-1;
    case 6
        totDATA=log2(log2(totDATA+1)+1);
        totDATA(totDATA>28)=28; %ct max
        if DATA.min_totDATA<0
            totDATA=totDATA-DATA.min_totDATA; % shift to ct min = 0
        end
        log2Ex=totDATA;
        DATA.max_totDATA=repmat(max(log2Ex),size(log2Ex,1),1);
        exponent=DATA.log_max_mRNA*(log2Ex./DATA.max_totDATA);
        totDATA=round(2.^exponent)-1;
        
end
end