

%% ANALYSE_EYE
%% Written by naveen at JLG on 12/10/18


clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');

cd(data_dir)

% LOADING FILE------------------------------------------------------- %% 
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE CELL DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile);
cd(PathName1);


i=10; FACTOR = 10000;

% EYE.EYEX = eye.X;
% EYE.EYEY = eye.Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of test %%%%%%%%%%%%%%%%%%%%%%%

EYE_Xvel = abs(diff(EYE.EYE_X{i,1}));
EYE_Yvel = abs(diff(EYE.EYE_Y{i,1}));

[p,pp] = findpeaks(EYE_Xvel,'Threshold',1);
EYE_Xb = EYE.EYE_X{i,1}; EYE_Yb = EYE.EYE_Y{i,1};

if length(pp)>=2
    
    if mod(length(pp),2)==0
        clear blink_pairs
        blink_pairs(1:length(pp)/2,1) = pp(1:2:length(pp));
        blink_pairs(1:length(pp)/2,2) = pp(2:2:length(pp));
        blink_pairs(1:length(pp)/2,3) = (blink_pairs(:,2)-blink_pairs(:,1));
        
        %%% removing blinks
        for j=1:length(pp)/2
            try
                EYE_Xb(blink_pairs(j,1)-20:blink_pairs(j,2)+20)=NaN;
                EYE_Yb(blink_pairs(j,1)-20:blink_pairs(j,2)+20)=NaN;
            catch
                EYE_Xb(blink_pairs(j,1):blink_pairs(j,2)+20)=NaN;
                EYE_Yb(blink_pairs(j,1):blink_pairs(j,2)+20)=NaN;
            end
        end
        %%% interpolating the removed values
        EYE_Xb = interpolate_NaN_n(EYE_Xb);
        EYE_Yb = interpolate_NaN_n(EYE_Yb);
        
        %%% smoothening
        EYE.EYE_Xn{i,1} = smooth(EYE_Xb,0.01,'rloess');
        EYE.EYE_Yn{i,1} = smooth(EYE_Yb,0.01,'rloess');
    else
        EYE.EYE_Xn{i,1} = [];
        EYE.EYE_Yn{i,1} = [];
    end
else
    EYE.EYE_Xn{i,1} = EYE.EYE_X{i,1};
    EYE.EYE_Yn{i,1} = EYE.EYE_Y{i,1};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of test %%%%%%%%%%%%%%%%%%%%%%%    
    
   color = [0 0 0];
   
if Infos(i,10)==1 colour = [0 0 1]; end
if Infos(i,10)==2 colour = [1 0 0]; end
    

EYE.EYE_X = EYE_POP_S{52, 1}.EYE_X;
EYE.EYE_Xn = EYE_POP_S{52, 1}.EYE_Xn;
EYE.EYE_Y = EYE_POP_S{52, 1}.EYE_Y;
EYE.EYE_Yn = EYE_POP_S{52, 1}.EYE_Yn;
i=14;

F = figure()
subplot(3,3,1)
hold on;
plot(EYE.EYE_X{i,1})
plot(EYE.EYE_Xn{i,1})
% ylim([-5 5])

subplot(3,3,2)
hold on;
plot(EYE.EYE_Y{i,1})
plot(EYE.EYE_Yn{i,1})
% ylim([-5 5])

subplot(3,3,3)
hold on;
plot(EYE.EYE_X{i,1}, EYE.EYE_Y{i,1})
plot(EYE.EYE_Xn{i,1}, EYE.EYE_Yn{i,1})
% xlim([-6-5 6+5]); ylim([-4-5 8+5])

subplot(3,3,4)
plot(EYE.EYE_Xn{i,1})
% ylim([-5 5])

subplot(3,3,5)
plot(EYE.EYE_Yn{i,1})
% ylim([-5 5])

subplot(3,3,6)
plot(EYE.EYE_Xn{i,1}, EYE.EYE_Yn{i,1})
% xlim([-6-5 6+5]); ylim([-4-5 8+5])




subplot(3,3,7)
plot(-200:200,EYE.EYE_Xn{i,1}(Infos(i,11)-200:Infos(i,11)+200),'color',colour)
% ylim([-5 5])

subplot(3,3,8)
plot(-200:200,EYE.EYE_Yn{i,1}(Infos(i,11)-200:Infos(i,11)+200),'color',colour)
% ylim([-5 5])

subplot(3,3,9)
plot(EYE.EYE_Xn{i,1}(Infos(i,11)-200:Infos(i,11)+200), EYE.EYE_Yn{i,1}(Infos(i,11)-200:Infos(i,11)+200),'color',colour)
% xlim([-6-5 6+5]); ylim([-4-5 8+5])




cd(Results_dir)
filename = strcat(FileName1,'-Trial-',num2str(i));
print(F, '-dpdf', filename, '-r400')







% % 
% % F = figure()
% % subplot(2,2,1)
% % hold on;
% % N = find(Infos(:,10)==1);
% % for i=1:25
% %     ii = N(i);
% %     if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
% %         plot(EYE.EYE_Xn{ii,1}, EYE.EYE_Yn{ii,1});
% %     end
% % end
% % xlim([-6-5 6+5]); ylim([-4-5 8+5])
% % grid on; grid minor
% % 
% % 
% % subplot(2,2,2)
% % hold on;
% % clear N
% % N = find(Infos(:,10)==2);
% % %  N = [N(3:7)];
% % for i=1:length(N)
% %     ii = N(i);
% %     if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
% %         plot(EYE.EYE_Xn{ii,1}, EYE.EYE_Yn{ii,1});
% %     end
% % end
% % xlim([-6-5 6+5]); ylim([-4-5 8+5])
% % grid on; grid minor;
% % ax.Xtick=[];
% % 
% % 




%EPOCH = Infos(ii,11)-200:Infos(ii,11)+600;


F = figure()

subplot(2,2,1)
hold on;
N = find(Infos(1:CHANGE,10)==1);
for i=1:45
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        %     if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
        plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
        %     end
    end
end
xlim([-6-5 6+5]); ylim([-4-5 8+5])
% grid on; grid minor
box on;

clear N
N = find(Infos(1:CHANGE,10)==2);
for i=1:5
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
            plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'r');
        end
    end
end




subplot(2,2,2)
hold on;
N = find(Infos(CHANGE+1:CHANGE+50,10)==1)+CHANGE;
for i=1:length(N)
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
    end
end
xlim([-6-5 6+5]); ylim([-4-5 8+5])
% grid on; grid minor
box on;

clear N
N = find(Infos(CHANGE+1:CHANGE+50,10)==2)+CHANGE;
for i=1:length(N)
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
            plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'r');
    end
end


COLOR_B = [25 188 157]/255; 
COLOR_R = [255 85 207]/255; 

subplot(2,2,3)
hold on;
N = find(Infos(1:CHANGE,10)==1);
for i=1:45
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        plot(-200:200,EYE.EYE_Xn{ii,1}(EPOCH),'-b')
        plot(-200:200,EYE.EYE_Yn{ii,1}(EPOCH),'color',COLOR_B)
%         plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
    end
end
ylim([-5 5]); xlim([-200 200]);


N = find(Infos(1:CHANGE,10)==2);
for i=1:5
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        plot(-200:200,EYE.EYE_Xn{ii,1}(EPOCH),'-r')
        plot(-200:200,EYE.EYE_Yn{ii,1}(EPOCH),'color',COLOR_R)
    end
end
ylim([-5 5]); xlim([-200 200]);





subplot(2,2,4)
hold on;
N = find(Infos(CHANGE+1:CHANGE+50,10)==1)+CHANGE;
for i=1:length(N)
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        plot(-200:200,EYE.EYE_Xn{ii,1}(EPOCH),'-b')
        plot(-200:200,EYE.EYE_Yn{ii,1}(EPOCH),'color',COLOR_B)
%         plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
    end
end
ylim([-5 5]); xlim([-200 200]);


N = find(Infos(CHANGE+1:CHANGE+50,10)==2)+CHANGE;
for  i=1:length(N)
    ii = N(i);
    EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
    if ~isempty(EYE.EYE_Xn{ii,1})
        plot(-200:200,EYE.EYE_Xn{ii,1}(EPOCH),'-r')
        plot(-200:200,EYE.EYE_Yn{ii,1}(EPOCH),'color',COLOR_R)
    end
end
ylim([-5 5]); xlim([-200 200]);



cd(Results_dir)
filename = strcat(NOME,'-OT_N');
print(F, '-dpdf', filename, '-r400')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % 
% % % 
% % % F = figure()
% % % subplot(2,2,1)
% % % hold on;
% % % EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % % eyex = NaN(length(N),length(EPOCH)); eyey=NaN(length(N),length(EPOCH));
% % % N = find(Infos(:,10)==1);
% % % for i=1:35
% % %     ii = N(i);
% % %     EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % %     if ~isempty(EYE.EYE_Xn{ii,1})
% % %         %     if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
% % %         eyex(ii,:) =  EYE.EYE_Xn{ii,1}(EPOCH);
% % %         eyey(ii,:) =  EYE.EYE_Yn{ii,1}(EPOCH);
% % % 
% % % %         plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
% % %         %     end
% % %     end
% % % end
% % % 
% % % 
% % % eyex_all(1,:) = nanmean(eyex);
% % % eyex_all(2,:) = nanstd(eyex);
% % % 
% % % eyey_all(1,:) = nanmean(eyey);
% % % eyey_all(2,:) = nanstd(eyey);
% % % 
% % % F = figure();
% % % subplot(2,2,1)
% % % errorbar(eyex_all(1,:),eyex_all(2,:))
% % % ylim([-5 5])
% % % 
% % % subplot(2,2,3)
% % % errorbar(eyey_all(1,:),eyey_all(2,:))
% % % ylim([-5 5])
% % % 
% % % 
% % % subplot(2,2,2)
% % % plot(eyex_all(1,:),eyey_all(1,:))
% % % ylim([-5 5]); xlim([-5 5])
% % % 
% % % 
% % % 
% % % 
% % % xlim([-6-5 6+5]); ylim([-4-5 8+5])
% % % grid on; grid minor
% % % 
% % % 
% % % % subplot(2,2,1)
% % % EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % % eyex = NaN(length(N),length(EPOCH)); eyey=NaN(length(N),length(EPOCH));
% % % N = find(Infos(:,10)==2);
% % % for i=1:length(N)
% % %     ii = N(i);
% % %     EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % %     if ~isempty(EYE.EYE_Xn{ii,1})
% % %         %     if length(EYE.EYE_Xn{ii,1})==length(EYE.EYE_Yn{ii,1})
% % %         eyex(ii,:) =  EYE.EYE_Xn{ii,1}(EPOCH);
% % %         eyey(ii,:) =  EYE.EYE_Yn{ii,1}(EPOCH);
% % % 
% % % %         plot(EYE.EYE_Xn{ii,1}(EPOCH), EYE.EYE_Yn{ii,1}(EPOCH),'-b');
% % %         %     end
% % %     end
% % % end
% % % 
% % % 
% % % eyex_all(1,:) = nanmean(eyex);
% % % eyex_all(2,:) = nanstd(eyex);
% % % 
% % % eyey_all(1,:) = nanmean(eyey);
% % % eyey_all(2,:) = nanstd(eyey);
% % % 
% % % F = figure();
% % % subplot(2,2,1)
% % % errorbar(eyex_all(1,:),eyex_all(2,:),'-r')
% % % ylim([-5 5])
% % % 
% % % subplot(2,2,3)
% % % errorbar(eyey_all(1,:),eyey_all(2,:),'-r')
% % % ylim([-5 5])
% % % 
% % % 
% % % subplot(2,2,2)
% % % plot(eyex_all(1,:),eyey_all(1,:),'-r')
% % % ylim([-5 5]); xlim([-5 5])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % 
% % % % %% HEATMAP GRID
% % % % 
% % % % Xgrid = -5:0.1:5;
% % % % Ygrid = -5:0.1:5;
% % % % 
% % % % nround = 2;
% % % % 
% % % % [X,Y] = meshgrid(Xgrid,Ygrid);
% % % % X = round(X,nround); Y = round(Y,nround);
% % % % X = X(:); Y = Y(:);
% % % % 
% % % % 
% % % % F = figure();
% % % % 
% % % % twodplotval = zeros(length(Xgrid)*length(Ygrid),1);
% % % % N = find(Infos(1:CHANGE,10)==1);
% % % % for i=1:length(N)
% % % %     ii = N(i);
% % % %     if ~isempty(EYE.EYE_Xn{ii,1})
% % % %         EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % % %         SIGNAL_x = EYE.EYE_Xn{ii,1}(EPOCH);
% % % %         SIGNAL_y = EYE.EYE_Yn{ii,1}(EPOCH);
% % % %         
% % % %         count=0;
% % % %         for n=1:length(SIGNAL_x)
% % % %             xval = double(round(SIGNAL_x(n),nround));   yval = double(round(SIGNAL_y(n),nround));
% % % %             ind = intersect(find(abs(X - xval) < 10^-3),find(abs(Y - yval) < 10^-3));
% % % %             twodplotval(ind) = twodplotval(ind)+1;
% % % %             count=count+1;
% % % %         end
% % % %         
% % % %         twodplotval = reshape(twodplotval,length(Xgrid),length(Ygrid));
% % % %     end
% % % %     
% % % %     subplot(2,2,2)
% % % %     hold on;
% % % %     plot(SIGNAL_x,SIGNAL_y);
% % % %     xlim([-5 5]); ylim([-5 5]);
% % % %     set(gca,'YDir','reverse');
% % % %     
% % % % end
% % % % 
% % % % twodplotval = twodplotval/length(N);
% % % % subplot(2,2,1)
% % % % imagesc(twodplotval/length(N));
% % % % colormap hot
% % % % 
% % % % 
% % % % %%%%%%%%%%%
% % % % 
% % % % twodplotval = zeros(length(Xgrid)*length(Ygrid),1);
% % % % N = find(Infos(1:CHANGE,10)==2);
% % % % for i=1:length(N)
% % % %     ii = N(i);
% % % %     if ~isempty(EYE.EYE_Xn{ii,1})
% % % %         EPOCH = Infos(ii,11)-200:Infos(ii,11)+200;
% % % %         SIGNAL_x = EYE.EYE_Xn{ii,1}(EPOCH);
% % % %         SIGNAL_y = EYE.EYE_Yn{ii,1}(EPOCH);
% % % %         
% % % %         count=0;
% % % %         for n=1:length(SIGNAL_x)
% % % %             xval = double(round(SIGNAL_x(n),nround));   yval = double(round(SIGNAL_y(n),nround));
% % % %             ind = intersect(find(abs(X - xval) < 10^-3),find(abs(Y - yval) < 10^-3));
% % % %             twodplotval(ind) = twodplotval(ind)+1;
% % % %             count=count+1;
% % % %         end
% % % %         
% % % %         twodplotval = reshape(twodplotval,length(Xgrid),length(Ygrid));
% % % %     end
% % % %     
% % % %     subplot(2,2,4)
% % % %     hold on;
% % % %     plot(SIGNAL_x,SIGNAL_y);
% % % %     xlim([-5 5]); ylim([-5 5]);
% % % %     set(gca,'YDir','reverse');
% % % %     
% % % % end
% % % % 
% % % % twodplotval = twodplotval/length(N);
% % % % subplot(2,2,3)
% % % % imagesc(twodplotval/length(N));
% % % % colormap hot
% % % % 
% % % % 
% % % % cd(Results_dir)
% % % % filename = 'EYE';
% % % % print(F, '-dpdf', filename, '-r400')
% % % % 
% % % % 
