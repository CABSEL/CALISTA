function [k_on, k_off, k_t] = get_rateconstants( k_tuples_mat,set_of_k )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% set of rate constants
k_on=set_of_k.on;
s_on=size(k_on);
k_on=reshape(k_on,1,[]);
k_off=set_of_k.off;
s_off=size(k_off);
k_off=reshape(k_off,1,[]);
k_t=set_of_k.t;
s_t=size(k_t);
k_t=reshape(k_t,1,[]);
k_d=set_of_k.d;
s_d=size(k_d);
k_d=reshape(k_d,1,[]);
nt=set_of_k.cell;


%number of cluster
noc=length(k_on);
% kon
k_tup_on=repmat(k_tuples_mat(:,1),1,noc);
diff=(k_tup_on-repmat(k_on,size(k_tup_on,1),1)).^2;
[~,in]=min(diff);
k_on=k_tup_on(in);
% koff
noc=length(k_off);
k_tup_off=repmat(k_tuples_mat(:,2),1,noc);
diff=(k_tup_off-repmat(k_off,size(k_tup_off,1),1)).^2;
[~,in]=min(diff);
k_off=k_tup_off(in);
% kt
noc=length(k_t);
k_tup_t=repmat(k_tuples_mat(:,3),1,noc);
diff=(k_tup_t-repmat(k_t,size(k_tup_t,1),1)).^2;
[~,in]=min(diff);
k_t=k_tup_t(in);

k_on=reshape(k_on,s_on(1),s_on(2));
k_off=reshape(k_off,s_off(1),s_off(2));
k_t=reshape(k_t,s_t(1),s_t(2));


end

