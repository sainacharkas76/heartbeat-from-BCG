function  [HR, pos_picchi, amp_picchi,Template] = SCG_template_matching_corr(axis,fs,tempAxis,template)

if isempty(template)
   
    Zbfilt_REV=axis;


    Zb30sec = {Zbfilt_REV};

    clc
    Template = []; 

    range_pre = 0.20;  

    range_post = 0.20;  


    blocchi = 1;

    for i=1:blocchi

            Zb10sec_part = Zb30sec{i,1};


            [massimo,indice] = max(Zb10sec_part(round(fs*2):round(fs*8),1)); % Finestra da 2 secondi agli 8 secondi. In questo modo si evita la presenza

            indice = indice + round(fs*2)-1;                      % di complessi non completi

            Template = cat(1,Template,Zb10sec_part(indice-round(range_pre*fs):indice+round(range_post*fs))); % Template di larghezza 400ms.

    end
else
    Template = template;
end

%% cross correlation calculation
clc


corr_Zbelly = {}; 


for i=1:blocchi

    Zb30sec_part = Zb30sec{i,1};

    corr_Zbelly_part = xcorr(Template(i,:),Zb30sec_part(:,1));

    corr_Zbelly_part = corr_Zbelly_part(1:length(Zb30sec_part),:);

    corr_Zbelly = cat(1,corr_Zbelly,corr_Zbelly_part);

end

%% thresholding Cross-Correlation
clc


Soglia = []; 


for i=1:blocchi

    Soglia = cat(1,Soglia,sqrt(mean((corr_Zbelly{i,1}).^2))*1.5); % Calcolo automatico della soglia (RMS*1.5)

    for j=1:length(corr_Zbelly{i,1})

        if (corr_Zbelly{i,1}(j,1) < Soglia(i,1))

            corr_Zbelly{i,1}(j,1) = 0;

        end

    end

end

%% refinement of the maximum values

clc

par_ricerca_cc = 40; 

for r = 1:blocchi

    for i= 1:length(round(corr_Zbelly{r,1}))

        if (corr_Zbelly{r,1}(i,1) > 0)

            pos = i;

            if (length(corr_Zbelly{r,1}) < pos+par_ricerca_cc)

                [massimo,indice] = max(corr_Zbelly{r,1}(pos:end,1));

                for j= pos:length(corr_Zbelly{r,1})

                    if (corr_Zbelly{r,1}(j,1) ~= massimo)

                        corr_Zbelly{r,1}(j,1)=0;

                    end

                end

                i=i+par_ricerca_cc;

            else

                [massimo,indice] = max(corr_Zbelly{r,1}(pos:pos+par_ricerca_cc,1));

                for j= pos:pos+par_ricerca_cc

                    if (corr_Zbelly{r,1}(j,1) ~= massimo)

                        corr_Zbelly{r,1}(j,1)=0;

                    end

                end

                i=i+par_ricerca_cc;

            end

        end

    end

end

%%
clc
% Elimination of the first value other than zero

for r= 1:blocchi

    for i= 1:length(corr_Zbelly{r,1})

        if (corr_Zbelly{r,1}(i,1) > 0)

            corr_Zbelly{r,1}(i,1) = 0;

            break;

        end

    end

end

figure
plot(corr_Zbelly{1, 1}  ) 
title('Z corr')

%% union of the signal and cross correlation

clc

Zbfilt_REV = []; 

corr_Zbfilt_REV = [];  

for i= 1:blocchi

    Zbfilt_REV = cat(2,Zbfilt_REV,Zb30sec{i,1});

    corr_Zbfilt_REV = cat(2,corr_Zbfilt_REV,corr_Zbelly{i,1});

end

%% better correlation

% resetting all negative values of the cross correlation

clc


for i=1:length(corr_Zbfilt_REV)

        if(corr_Zbfilt_REV(i,1) <= 0)

            corr_Zbfilt_REV(i,1) = 0;

        end

end


close all

figure(1)
plot (Template)
savefig(strcat('Template',tempAxis,'.fig'))
close all

%%

corr_Zbfilt_dritto = corr_Zbfilt_REV(end:-1:1,1);
picco = find(corr_Zbfilt_dritto~=0);
picco

% time threshold
picco_temp = picco;
count = 0;

for i = 1:size(picco_temp)-1
    if strcmp (picco_temp(i), NaN) ==1
        continue
    end
    if (picco_temp(i+1,1) - picco_temp(i,1)) < 0.4*fs
        picco_temp (i+1) = NaN;
    end
end
picco_temp = picco_temp(isnan(picco_temp)==0);

figure
plot (axis,'k')
hold on
plot (picco_temp, axis(picco_temp),'r*')

figure
plot (axis,'k')
hold on
plot (picco_temp, axis(picco_temp),'r*')
savefig(strcat('Picchi_TempMatch',tempAxis,'.fig'))
close all

pos_picchi = picco_temp;
amp_picchi = axis(picco_temp);


%% Compute HR

clc
HR = 60./(diff(picco_temp))*fs;

end

