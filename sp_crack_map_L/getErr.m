function out=getErr(a,n)
% Joshua Yeh
% Created: 18/5/29
% This function computes the error based on the number significant digits.
%% INPUT VARIABLES
% a: is the number rounded to n significant digits
% n: # of sig fig a was rounded to
% 
%% OUTPUT VARIABLES
% out: the error associated with a based on the number of significant
% digits
% 
txt=num2str(a,'%10.10f');

% find the decimal place
dec=[];
for dum=1:numel(txt)
    if strcmp(txt(dum),'.')
        dec=dum;
    end
end

% find significant digits
flag=1;
flag1=1;
sf=[];
count=0;
c=1;
while flag==1
% for dum=1:numel(txt)
    if count==n
        break
    end
    
    if c<=numel(txt)
        if (~strcmp(txt(c),'0')&&~strcmp(txt(c),'.'))&&flag1==1
            flag1=0;
            sf=cat(1,sf,c);
            count=count+1;
        elseif isstrprop(txt(c),'digit')&&flag1==0
            sf=cat(1,sf,c);
            count=count+1;
        end
        
    else
        sf=cat(1,sf,c);
        count=count+1;
    end
    c=c+1;
end

% decimal placement of least significant digit
placement=dec-sf(end);

if placement>0
    placement=placement-1;
end

out=5*10^(placement-1);
