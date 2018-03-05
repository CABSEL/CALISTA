function my_diff=get_diff(v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

len=length(v);
my_diff=zeros(1,len-1);
for i=1:len-1
    my_diff(i)=v(i+1)-v(i);

end

