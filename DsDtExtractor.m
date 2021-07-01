

for ii=1:80
    


[x,y]=find(sat_change_vs_time==max(sat_change_vs_time(:,ii)));
CellIndex(ii)=x;
[I,J] = ind2sub([1,9],x)
distance(ii)=0.05*((I-1)^2+(J-9)^2)^0.5
maxSw(ii)=max(sat_change_vs_time(:,ii));
meanSw(ii)=mean(sat_change_vs_time(:,ii));





end

OutPut(:,1)=distance';
OutPut(:,2)=maxSw';
plot(distance,maxSw,'o')
figure
plot(distance,meanSw,'*')



