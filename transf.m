function M=transf(M,method, dim,parms)
%     M=transf(M,method, dim,parms)
%Supported transformations:   [Z]
% Z zscore
% C center/de-mean
% F Fisher
% R Range (parms=[min max]) [0 1]
% B Binning  (parms=+/-nbin, negative calls quantile) [-10]
% A Asin(sqrt) 
%
% S extremevalues Suppress (atan?)
% G extremevalues Gainer (fisher?)
%3D support?
%ok for F A
%ok for Z C R B and dim=1


if ~nargin
    help transf
else

    
    


if ~exist('method','var')||isempty(method)
method='Z';
end
method=upper(method(1));

switch method
    %these are direct, population independent operations, no reshaping is needed
        case 'F'
        M(M>=1)=1-eps;    M(M<=-1)=-1+eps;           
        M=0.5*log((1+M)./(1-M));       
%      M=atanh(M);        
        case 'A'
        M=asin(sqrt(M));        
otherwise
    
    
if ~exist('dim','var') || isempty(dim)
    resh=true;
    flip=false;   
    [d1,d2]=size(M);
    M=M(:);      
elseif dim==2
    resh=false;    
    flip=true;    
    M=M';
    [d1,d2]=size(M);
else
    resh=false;    
    flip=false;
   [d1,d2,d3]=size(M);   
end     
    

switch method
case 'R'
        if ~exist('parms','var')
        extra=false;
        elseif all(parms==[0 1])
        extra=false;
        else
        extra=true;                
        rangev=range(parms); 
        minv=min(parms);
        end
case 'S' | 'G'
    method='R';  
    minv=-1; rangev=2;
case 'B'
         if ~exist('parms','var')
         b=quantile(M(:),linspace(0,1,10+1));
         b([1 end])=[-Inf Inf];
         else
         if numel(parms)==1
                if resh
                if parms>0
                     b=linspace(min(M(:)),max(M(:)),parms+1);
                else
                     b=quantile(M(:),linspace(0,1,-parms+1));                                             
                end
               b([1 end])=[-Inf Inf];                          
                else
                b=zeros(d2,abs(parms)+1);
                for n=1:d2
                if parms>0
                     b(n,:)=linspace(min(M(:,n)),max(M(:,n)),parms+1);
                else
                     b(n,:)=quantile(M(:,n),linspace(0,1,abs(parms)+1));                                             
                end                                       
                end
                 b(:,1)=-Inf;  b(:,end)=Inf; 
                end
         else
                b=parms;
         end
         end
end
    
   


if resh
switch method
    case 'Z'
        M=(M-nanmean(M))/nanstd(M);    
    case 'C'
        M=M-nanmean(M);        
    case 'R'
        M=M-min(M);
        M=M./max(M);
        if extra
        M=(M*rangev)+minv;
        end
    case 'B'
        [~,M]=histc(M,b);
end
M=reshape(M,d1,d2);
else
   switch method
    case 'Z'
%         M=zscore(M);
        M=(M-repmat(nanmean(M),[d1 1 1]))./repmat(nanstd(M),[d1 1 1]);
    case 'C'
        M=M-repmat(nanmean(M),[d1 1 1]);
    case 'R'        
        M=M-repmat(min(M),[d1 1 1]);
        M=M./repmat(max(M),[d1 1 1]);
        if extra
        M=(M*rangev)+minv;
        end
        
       case 'B'
        for m=1:d3   
        for n=1:d2
        [~,M(:,n,m)]=histc(M(:,n,m),b(n,:));
        end
        end
   end
end 
     

if flip
    M=M';
end
end
end