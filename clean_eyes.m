%% function to 1) remove blinks 2) remove artifacts from eye signals
% written by naveen at JLG on 9/13/19

function [EYE_X,EYE_Y] = clean_eyes(eyeX,eyeY)


LEN = size(eyeX,1);

%%%%%% remove random artifacts

% FACTOR = 10000;
FACTOR = 1;
Thr = 0.5;

parfor i=1:LEN
    i
    if ~isempty(eyeX{i,1})
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of test %%%%%%%%%%%%%%%%%%%%%%%
        
        EYE_Xvel = abs(diff(eyeX{i,1}));
        EYE_Yvel = abs(diff(eyeY{i,1}));
        
        [p,pp] = findpeaks(EYE_Xvel,'Threshold',Thr*FACTOR);
        EYE_Xb = eyeX{i,1}; EYE_Yb = eyeY{i,1};
        
        if length(pp)>=2
            %         if mod(length(pp),2)==0
            EYE_Xb(find(EYE_Xb>Thr*FACTOR)) = NaN;   EYE_Yb(find(EYE_Yb>Thr*FACTOR)) = NaN;
            try
                EYE_Xb = interpolate_NaN_n(EYE_Xb);
                EYE_Yb = interpolate_NaN_n(EYE_Yb);
                
                EYE_X{i,1} = smooth(EYE_Xb,0.01,'rloess');
                EYE_Y{i,1} = smooth(EYE_Yb,0.01,'rloess');
            catch
                EYE_X{i,1} = [];
                EYE_Y{i,1} = [];
            end
            
        else
            EYE_X{i,1} = eyeX{i,1};
            EYE_Y{i,1} = eyeY{i,1};
        end
    else
        EYE_X{i,1} = [];
        EYE_Y{i,1} = [];
    end
    
end





end