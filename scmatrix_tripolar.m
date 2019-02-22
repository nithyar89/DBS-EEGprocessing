function scmatrix(XY,M,sc,im,x,ev,i,lbl,vref)
%   scmatrix(    XY  ,M,sc,im,x,ev,i,lbl,vref)
%   scmatrix(chanlocs,M,sc,im,x,ev,i,lbl,vref)
if ~nargin
    help scmatrix
    return
end

% set(gcf,'color','k')
clf
%create figure layout
if isstruct(XY) %assuming XY is chanlocs/plocs
  try  
    [X,Y]=ch2plocs(XY); %chanlocs
    nn={XY.labels}';
  catch
     X=XY.x; Y=XY.y; %plocs
     nn=XY.lbl;
  end
X=transf(X,'R',1,[0 1]);
Y=transf(Y,'R',1,[0 1]);       
n=numel(X);
else
X=transf(XY(:,1),'R',1,[0 1]);
Y=transf(XY(:,2),'R',1,[0 1]);   
n=numel(X);   
nn=cell(n,1);
for j=1:n
    nn{j}=num2str(j);
end
end

if ~exist('i','var')||isempty(i) %subset of electrodes, ordinal
    i=1:n;
else
    i2=false(n,1); i2(i)=true; i=find(i2)';
end

if exist('lbl','var') && ~isempty(lbl)
    lblon=1;
else
    lblon=0;
end

if n>2
w=min(pdist([X(i(:)) Y(i(:))]))/1.4143; %for square inlays    
X=transf(X,'R',1,[w/2 1-w/2]);
Y=transf(Y,'R',1,[w/1.5 1-w/1.5]);   
w=1.2*min(pdist([X(i(:)) Y(i(:))]))/sqrt(2); %for square inlays
elseif n==2
w=.5;  
X=transf(X,'R',1,[w/2 1-w/2]);
Y=transf(Y,'R',1,[w/2 1-w/2]);   
w=min(pdist([X Y]))/sqrt(2); %for square inlays
elseif n==1
    X=.5; Y=.5; w=0.5;
end
p.npf=n;     
p.maximsperfig=n;
p.lm=X-(w/2);  p.bm=Y-(w/2); %shifted
p.w=repmat(w,[n 1]);
p.h=repmat(w,[n 1]);

   


%handle input
dim=size(M);
eldim=find(dim==n,1,'first'); %the first if ambiguous
if ~exist('im','var') %imagesc/plot choice is not defined
if numel(dim)==3
    im=true; %it will be imagesc
else
    im=false; %it will be a plot
end
end

if ~exist('sc','var')||isempty(sc)
     sc=quantile(M(~isnan(M)),[.05  .95]);
     if all(~sc)
         sc=[-1 1];
     end
elseif numel(sc)==1
    sc=[-abs(sc) abs(sc)];
end

if ~exist('x','var')||isempty(x)
    dim2=dim; dim2(eldim)=[];
     x=1:max(dim2);
end

showev=false;
if exist('ev','var') && ~isempty(ev)
     showev=true;
end




hh=zeros(n,1);
if im    
for j=i
hh(n)=subplot('Position',[p.lm(j) p.bm(j) p.w(j) p.h(j)]);
switch eldim
    case 3
    IM=M(:,:,j);
    case 2
    IM=squeeze(M(:,j,:));        
    case 1
    IM=squeeze(M(j,:,:));                
end
if size(IM,2)~=numel(x)
    IM=IM';
end
yt=size(IM,1);
imagesc(x,1:yt,IM,sc)
axis xy
text(median(x),yt*1.1,nn{j},'VerticalAlignment','top','HorizontalAlignment','center','color',[.4 .4 .9])    

if showev
        line(repmat(ev(:),[1 2])',[0 yt+2],'Color',[0 0 0])
        text(ev(:),repmat(sc(1),[1 numel(ev)]),num2str(ev(:)),'VerticalAlignment','top','HorizontalAlignment','center','color',[1 1 1]*.6)    
end
if exist('vref','var')
text(ones(numel(vref),1)*x(1)-70,vref(:),num2str(vref(:)),'VerticalAlignment','top','HorizontalAlignment','center','color',[1 1 1])    
end
axis off
end

else
%plot
yl=range(sc)*.1;
yl=[sc(1)-yl sc(2)+yl];
xl=range(x)*.1;
xl=[x(1)-xl x(end)+xl];
[~,center_el]=min(pdist2(mean([X Y]),[X Y]));

for j=i
hh(n)=subplot('Position',[p.lm(j) p.bm(j) p.w(j) p.h(j)]);

axis([xl yl])    
switch eldim
    case 3
    y=M(:,:,j);
    case 2
    y=squeeze(M(:,j,:));        
    case 1
    y=squeeze(M(j,:,:));                
end
if j==1 && lblon
    plot(x,y,'Linewidth',2)
    legend(lbl)
    axis([xl yl])    
end
hold on
line(xl,[0 0],'Color',[.4 .4 .9]*1.1,'LineStyle','-')
line(xl,[sc(1) sc(1)],'Color',[1 1 1]*.6,'LineStyle',':')
line(xl,[sc(2) sc(2)],'Color',[1 1 1]*.6,'LineStyle',':')

text(double(median(x)),double(yl(2)),nn{j},'HorizontalAlignment','center','color',[.4 .4 .9])
plot(x,y,'Linewidth',2)

    if j==center_el
    text(double(xl(1)),double(sc(1)),num2str(sc(1)),'VerticalAlignment','Middle','HorizontalAlignment','Right','color',[1 1 1]*.6)    
    text(double(xl(1)),double(sc(2)),num2str(sc(2)),'VerticalAlignment','Middle','HorizontalAlignment','Right','color',[1 1 1]*.6)    
    end


if showev
%         line(repmat(ev(:),[1 2])',[sc(1) sc(2)],'Color',[.4 .4 .9])
        line(repmat(ev(:),[1 2])',[sc(1) sc(2)],'LineStyle',':','color',[.4 .4 .4])  
        if j==center_el
        text(ev(:),repmat(sc(1),[1 numel(ev)]),num2str(ev(:)),'VerticalAlignment','top','HorizontalAlignment','center','color',[1 1 1]*.6)    
        end
end
axis off
end
end







    
    


    
    





    


