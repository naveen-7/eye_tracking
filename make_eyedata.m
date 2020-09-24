
%% function to align the eye data to an event
% writtne by naveen at JLG on 3/26/20

function OUTPUT = make_eyedata(MAT,ALIGN,LIM1,LIM2)

LEN = length(LIM1:LIM2);
OUTPUT = NaN(length(ALIGN),LEN);
Len = cellfun(@length, MAT);

for i=1:length(ALIGN)
    try
        clear dum temp;
        dum = -ALIGN(i):Len(i)-ALIGN(i)-1;
        temp = MAT{i,1}(find(dum==LIM1):find(dum==LIM2));
        temp = temp-nanmean(temp);
        if abs(nanmean(temp(1:5)))<=10
            OUTPUT(i,:) = temp;
        else
            OUTPUT(i,:) = NaN;
        end
    end
end

% OUTPUT = nanmean(OUTPUT);

end