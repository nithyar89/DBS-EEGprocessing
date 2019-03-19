function M=struct2mat(S,field,all)
%           M=struct2mat(S,field,all)
%It converts directly structs into matrix by extracting the content of
%the specified 'field' from 'S'.
%If field is not specified all fields are returned. If fields have mixed
%data types only numerical data are kept. 
%'all' specify how to handle the case of field contents of different
%lengths: if 0, the shortest common length is used, otherwise shorter
%lengths are filled with NaNs [1]. 

if ~nargin
    help struct2mat
else
if ~exist('all','var')
    all=1;
end

if ~isstruct(S)
error('Input to struct2mat is not a structure!')
else
    
    if exist('field','var') && ~isempty(field)
        loc=strcmp(field,fieldnames(S));
        S=struct2cell(S);
        S=S(loc,:,:,:,:,:,:,:,:);        
    else
            S=struct2cell(S);
    end
    try
    M=squeeze(cell2mat(S));
    catch ME
     if strcmp(ME.identifier, 'MATLAB:cell2mat:UnsupportedCellContent')
            M=squeeze(S);
     elseif ischar(S{1}) && strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') 
            M=squeeze(S);
     else
        if strcmp(ME.identifier,'MATLAB:cell2mat:MixedDataTypes')
            nf=size(S,1);
            ft=zeros(1,nf);
            for n=1:nf
                ft(n)=ischar(cell2mat(S(n)));
            end
                S=S(~ft,:,:,:,:);        
        end
        
        
        
        S=S(:);
        l=zeros(size(S,1),1);
        for n=1:size(S,1)
            l(n)=numel(S{n});
        end
        if all
           l=max(l);
           M(1:l,1:size(S,1))=NaN;
           for n=1:size(S,1)
               c=cell2mat(S(n));   
               M(1:numel(c),n)=c;
           end
        else
           l=min(l);
           M=zeros(l,size(S,1));
           for n=1:size(S,1)
               c=cell2mat(S(n));   
               M(:,n)=c(1:l);
           end
        end 
     end
    end
    if isvector(M)
        M=M(:);
    end
end
end