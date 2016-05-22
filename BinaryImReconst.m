function [Ic]=BinaryImReconst(Im,options)
%Inputs:
%    Im: image with missing data
%    ws: window size
%    lambda: weight of model (vs data)
%    step: step determines the overlap (step=1=> full overlap)
%         1<step<ws
%Outputs:
%    Ic: completed image
step=options.step;
ws=options.windowsize;
lambda=options.lambda;
eta=options.eta;
eps=options.eps;
nmaxiter=options.gdmaxiters;
number_of_iters=options.iters;
% options.iters=4;

is=size(Im);%image size
N=prod(ws);%# of pixels in a window
Ic=zeros(is);

n1=ws(1);n2=ws(2);
x1=-(n1-1)/2:(n1-1)/2;
x2=-(n2-1)/2:(n2-1)/2;
[x1,x2]=meshgrid(x1,x2);
x1=x1(:)';
x2=x2(:)';
X=[x1;x2];


for i=1:number_of_iters%options.iters
    Icnew=zeros(is);
    Is=Icnew;
    for x=1:step:is(1)-ws(1)+1
        waitbar((i-1+x/(is(1)-ws(1)+1))/options.iters);
        for y=1:step:is(2)-ws(2)+1
            pm=Im(x:x+ws(1)-1,y:y+ws(2)-1);% the patch with missing data
            p=Ic(x:x+ws(1)-1,y:y+ws(2)-1);% the patch from completed image
            pmv=pm(:);
            ind=find(~isnan(pmv));
            q=pmv(ind);
            %% W-step:
            for j=1:1
                if strcmp(options.model,'ULEM')
                    if (i==1)&&(j==1)
                        [~,~,c]=UnconstSingleGD(X(:,ind),pm(ind)',eta,eps,nmaxiter,X);
                    else
                        [~,~,c]=UnconstSingleGD(X,p(:)',eta,eps,nmaxiter,X);
                    end
                elseif strcmp(options.model,'UCM')
                    if (i==1)&&(j==1)
                        [~,~,c]=UnconstCombGD(X(:,ind),pm(ind)',eta,eps,nmaxiter,X);
                    else
                        [~,~,c]=UnconstCombGD(X,p(:)',eta,eps,nmaxiter,X);
                    end
                elseif strcmp(options.model,'CCM')
                    if (i==1)&&(j==1)
                        [~,~,c]=ConstCombGD_fix_channel(X(:,ind),pm(ind)',eta,eps,nmaxiter,X);
                    else
                        [~,~,c]=ConstCombGD_fix_channel(X,p(:)',eta,eps,nmaxiter,X);
                    end
                end
                %% P-step:
                H=eye(N);
                H=H(ind,:);
                
                
                pnew=(H'*H+lambda*eye(N))\(H'*q+lambda*c');%%%%%%%%%%%%%%%%%cccccc
                p=pnew;
            end
            Icnew(x:x+ws(1)-1,y:y+ws(2)-1)=Icnew(x:x+ws(1)-1,y:y+ws(2)-1)+reshape(pnew,ws);
            Is(x:x+ws(1)-1,y:y+ws(2)-1)=Is(x:x+ws(1)-1,y:y+ws(2)-1)+ones(ws);
            %imshow(Icnew./Is);
            %pause
            
        end
    end
    Icnew=Icnew./Is;
%     Ic=double(Icnew>0.5);% in this part we dedicate threshold of finded pic for next iter
%     Ic=Icnew;
    Ic_threshold=double(Icnew>0.5);
     Ic=Ic_threshold;
%     subplot(1,options.iters,i);
%     figure
%     imshow(Ic)
    figure
    imshow(Ic_threshold)
    title('thresholded result image');
    %keyboard
%     
%     % calculating performance
%     load ti_channel
%     ti=rot90(ti,-2);
%     ti= flipdim(ti,2);
% 
%     Iorg=ti;
%     sub_Iorg_and_Ic_threshold=abs(Iorg-Ic_threshold);
%     sum_sub_Iorg_and_Ic_threshold=sum(sum(sub_Iorg_and_Ic_threshold));
%     [q40,q41]=size(ti);
%     performance(1,i)=100-(abs(sum_sub_Iorg_and_Ic_threshold)/(q40*q41))*100
  
    
end