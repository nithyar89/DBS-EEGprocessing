function     [out1,y,xyz,out4]=ch2plocs(chanlocs,in2,m)
%         plocs=ch2plocs(chanlocs)
%     [x,y,xyz]=ch2plocs(chanlocs)
%         plocs=ch2plocs(chanlocs,gridsc)
%         plocs=ch2plocs(chanlocs,AREA)
%     [x,y,xyz]=ch2plocs(chanlocs,gridsc)
%     [XX,YY,ZZ,MM]=ch2plocs(chanlocs,gridsc,m)


if ~nargin
    help ch2plocs
else
            Th=deg2rad(cell2mat({chanlocs.theta})');
            Rd=cell2mat({chanlocs.radius})';
            [x,y] = pol2cart(Th,Rd);
            nrad=pi/2; %rotate of nrad
            nrad=pi+pi/2;
            y=-y;
            allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*nrad);
            
            if ~exist('in2','var') 
            plocs.x = imag(allcoords);
            plocs.y = real(allcoords);          
            plocs.xyz=[struct2mat(chanlocs,'X')  struct2mat(chanlocs,'Y') struct2mat(chanlocs,'Z')];
            plocs.lbl={chanlocs.labels}';
            if nargout>1
               out1=plocs.x;                
               y=plocs.y;
               xyz=plocs.xyz;        
            else
                out1=plocs;
            end            
            else
                if isstruct(in2) %AREA
                plocs.x = el2a(imag(allcoords),in2,1);
                plocs.y = el2a(real(allcoords),in2,1);
                plocs.xyz=el2a([struct2mat(chanlocs,'X')  struct2mat(chanlocs,'Y') struct2mat(chanlocs,'Z')],in2,1);                            
                plocs.lbl={in2.name}';                
                if nargout>1
                   out1=plocs.x;                
                   y=plocs.y;
                   xyz=plocs.xyz;        
                else
                   out1=plocs;
                end                
                else
                if nargin==2    
                plocs.gridsc=in2;            
                plocs.x = imag(allcoords);
                plocs.y = real(allcoords);          
                plocs.xyz=[struct2mat(chanlocs,'X')  struct2mat(chanlocs,'Y') struct2mat(chanlocs,'Z')];                
                plocs.lbl={chanlocs.labels}';                
                X=nanmedian(surfint(plocs.x,chanlocs,plocs.gridsc));
                Y=nanmedian(surfint(plocs.y,chanlocs,plocs.gridsc),2);
                i=1:numel(X);
                plocs.X=interp1(i(~isnan(X)),X(~isnan(X)),i,'linear','extrap');
                plocs.Y=interp1(i(~isnan(Y)),Y(~isnan(Y)),i,'linear','extrap');    
                if nargout>1
                   out1=plocs.x;                
                   y=plocs.y;
                   xyz=plocs.xyz;        
                else
                   out1=plocs;
                end 
                else
                    
                    
                    
                    
                end
                end
            end
            
      

end        
