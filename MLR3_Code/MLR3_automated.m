%% 1. Loading necessary mat files, Raw Counts to TPM and RNA to microarray conversion
clear;
clc;

%cd  F:/IISC/EMT/EMT_Score/EMT_Score_Rcodes/EMT_Scoring_RNA_seq/MLR3_Code                                      %commented when using in ubuntu


load('RelevantData.mat'); 
dirInfo = importdata("matlab_input_temp.txt");
cd  ../Data_generated/
	

	MLR3data = readtable(char(dirInfo(1,1))); % Predictors+Normalizers
	data = table2array(MLR3data(1:25,2:end));

	idx = readtable(char(dirInfo(2,1)));
	NCI_data = DataNCI60(table2array(idx), :);
	NCI_data(:, 34) = [];
	new_norm = mean(mean(data(6:25,:)));
	NCI_norm = mean(mean(NCI_data(6:25,:)));

	d = new_norm - NCI_norm;
	Norm_data = data - d;

	a1 = 1; a2 = 2; a3 = 3;
	c1 = 1; c2 = 2; c3 = 3;

    
cd ../MLR3_Code/
    
%Predictions
%% 2. (CLDN7, VIM/CDH1)

	NUM = size(data, 2);

	Data = data;
	NormData = Norm_data;

	[YhatDataEMTGenes8pair1, PredictionsDataEMTGenes8pair1] = MLR3(...
		B1, GeneList1,...
		[{'CLDN7'}; {'VIM/CDH1'}],...
	[Data(a1,:); Data(a2,:)./Data(a3,:)]);

	[YhatDataEMTGenes8pair1Norm, PredictionsDataEMTGenes8pair1Norm] = MLR3(...
        	B1, GeneList1,...
        	[{'CLDN7'}; {'VIM/CDH1'}],...
        [NormData(a1,:); NormData(a2,:)./NormData(a3,:)]);

	YHat1 = YhatDataEMTGenes8pair1; YHat1_norm = YhatDataEMTGenes8pair1Norm;

	ScoreEMT1 = nan(NUM, 1);
		for j = 1:NUM
			if YHat1(j, 1) > YHat1(j, 3)
            		ScoreEMT1(j) = YHat1(j, 2);
        	else
            		ScoreEMT1(j) = 2.0 - YHat1(j, 2);
        	end
    	end

    	ScoreEMT1_norm = nan(NUM, 1);
    		for j = 1:NUM
        		if YHat1_norm(j, 1) > YHat1_norm(j, 3)
            		ScoreEMT1_norm(j) = YHat1_norm(j, 2);
        	else
            		ScoreEMT1_norm(j) = 2.0 - YHat1_norm(j, 2);
        	end
    	end

%Pre-normalization and Post-normalization scatterplots

cd './Plots/'

	gseID = split(string(dirInfo(1,1)),"_")
	imageName = gseID(1) + "_scoreEMT1.png"

	figure(1);
	plot(0,0,'.k'); plot(0,0,'.b'); plot(0,0, '.r');
	scatter(DataNCI60(c1,:),log2((DataNCI60(c2,:)./DataNCI60(c3,:)+1)),'k*');
	hold on
	scatter(data(a1,:),log2((data(a2,:)./data(a3,:)+1)),'b*');
	hold on
	scatter(Norm_data(a1,:),log2((Norm_data(a2,:)./Norm_data(a3,:)+1)),'r*');

	xlabel('CLDN7 Expression'); ylabel('VIM/CDH1 Expression');
	legend('NCI60','Data','NormData');
	title('CLDN7 vs VIM/CDH1');
	hold off

	saveas(gcf,imageName);



%% 5. Compiling and saving scores and correlations
cd  ../../Data_generated/

	Results = table(ScoreEMT1, ScoreEMT1_norm);
	writetable(Results,string(dirInfo(3,1)),'Delimiter','\t');
