function cmap=TwoColorSpanningCmap(color1,color2,cmapLength)
    if(nargin<2)
        cmapLength=128;
    end
   cmap=zeros(cmapLength,3);
   for i=1:3
    cmap(:,i)=linspace(color1(i),color2(i),cmapLength);
   end
end
