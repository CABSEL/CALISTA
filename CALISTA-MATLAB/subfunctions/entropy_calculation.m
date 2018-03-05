function [ H] = entropy_calculation( DDD,min_num_cells )

for i=1:size(DDD,2)
    [n_genes,n_cells]=size(DDD{i});
    for j=1:n_genes
        x = DDD{i}(j,:);
        sorted_x=sort(x);
        in_bin=[];
        for xx=1:length(sorted_x)
            if xx==1
                num_bin=1;
                in_bin(num_bin)=0;
            end
            in_bin(num_bin)=in_bin(num_bin)+1;
            
            if xx<length(sorted_x) & sorted_x(xx)~=sorted_x(xx+1) & in_bin(num_bin)>=min_num_cells
                num_bin=num_bin+1;
                in_bin(num_bin)=0;
            end
        
        end  
        binning=[];
        binning=mat2cell(sorted_x+1,1,in_bin);
        min_mRNA_in_bin=1;
        for equiBIN=1:length(in_bin)
            max_mRNA_in_bin=max(binning{equiBIN});
           
%             if max_mRNA_in_bin==0
%                 max_mRNA_in_bin=1;
%             end
            in_bin(2,equiBIN)=length(min_mRNA_in_bin:max_mRNA_in_bin);
             min_mRNA_in_bin= max_mRNA_in_bin+1;
        end
        
        frequency=in_bin(1,:)/n_cells;
        bin_size=in_bin(2,:);
        H(i,j) = -sum(frequency.*log2(frequency./bin_size));
            
    end


end
