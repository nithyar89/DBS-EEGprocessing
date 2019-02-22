function [contents] = get_directory_contents(path, pattern1, pattern2)
%Get Directory Contents - get list of strings corresponding to files in a
%given directory
%   path = variable containing string of directory
%   pattern1 = starting pattern in coded files (e.g. 'subj')
%   pattern2 = file extension (e.g. '.mat')
%   contents = matrix cell of directory files

% version 1.0, 14/11/2014

% copyright Nik Murphy, 2014, n.murphy@ncl.ac.uk

%get directory information
listing = dir(path);
%extract file names
si = size(listing,1);
s=cell(si,1);
for iter = 1:si
    s{iter}=listing(iter,1).name;
end
%find names with appropriate pattern
d = cell(si,1);
pattern = pattern1;
for iter = 1:si
    b = s{iter};
    c = strfind(b,pattern);
    if c >=1
        d{iter}=b;
    end
end
%remove inappropriate results
out = zeros(1,si);
for iter = 1:si
    t = isempty(d{iter});
    if t==1
        out(iter)=iter;
    end
end
for iter =1:si
    if out(iter)>0
        out(iter) = 1;
    end
end
out = logical(out);
d(out,:)=[];
out(:,:)=[];
%find names with correct file extension
si = size(d,1);
e = cell(si,1);
pattern = pattern2;
for iter = 1:si
    b = d{iter};
    c = strfind(b,pattern);
    if c >=1
        e{iter}=b;
    end
end
%remove inappropriate results
% out = zeros(1,si);
for iter = 1:si
    t = isempty(e{iter});
    if t==1
        out(iter)=iter;
    end
end
% test = isempty(out);
% if test == 0
%     [out]=removeZeroColumnsByElem(out); %see subroutine
%     e(out,:)=[];
% end
R=(~cellfun('isempty',e)); % was getting issues with above routine.  
% cell function works more efficiently.
%output
contents = e(R);
end
%% Subroutine

% function [ out ] = removeZeroColumnsByElem( matrix )
%
%
% inputs = matrix;
% realinputs = [];
% %get only the ones that actually have data
% count = 1;
% for i = 1:size(inputs,2)
%
%     % are all the values in this column zero?
%
%     isZero = 1;
%
%     for j = 1:size(inputs,1)
%         if inputs(j,i) ~= 0
%             isZero = 0;
%
%         end
%     end
%
%     if isZero == 0
%
%         realinputs(:,count) = inputs(:,i);
%         count = count+1;
%     end
%
% end
%
% out = realinputs;
%
% end
