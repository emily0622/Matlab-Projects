matrix_struct = load("PS01_dataSet/wordVecV.mat");
W = transpose(matrix_struct.V);
results_dis = zeros(10,10);
final_result_dis = inf;
final_result_indicies_dis = [0,0];
results_ang = zeros(10,10);
final_result_ang = 360;
final_result_indicies_ang = [0,0];

%  0 0 0 0
%  0 0 0 0
%  0 0 0 0
%  0 0 0 0 
% (2,1) 
% (3,1) (3,2)
% (4,1) (4,2) (4,3)

% 1 B. J. Cole
% 2 Mary J. Blige
% 3 Jessica Feshbach
% 4 Susie Au
% 5 Geoff Brown (baseball)
% 6 John Holland (composer)
% 7 James Forder
% 8 Public image of George W. Bush
% 9 Barack Obama
% 10 George W. Bush



for results_column = 2:10;
    for results_row = 1:(results_column-1);
        
        %part a part 1
        u = W(results_column,:);
        v = W(results_row,:);
        
        D = norm(u-v);
        results_dis(results_column,results_row) = D;
        
        if final_result_dis > D
            final_result_dis = D;
            final_result_indicies_dis = [results_column,results_row];
        end
        
        %part a part 2
        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));
        results_ang(results_column,results_row) = ThetaInDegrees;
        
        if final_result_ang > ThetaInDegrees
            final_result_ang = ThetaInDegrees;
            final_result_indicies_ang = [results_column,results_row];
        end
        
    end
end


results_dis
final_result_dis
final_result_indicies_dis

results_ang
final_result_ang
final_result_indicies_ang


for results_column = 2:10;
    for results_row = 1:(results_column-1);
        
        %part b part 1
        u_1 = W(results_column,:);
        u = u_1/(norm(u_1));
        v_1 = W(results_row,:);
        v = v_1/norm(v_1);
        
        
        D = norm(u-v);
        results_dis(results_column,results_row) = D;
        
        if final_result_dis > D
            final_result_dis = D;
            final_result_indicies_dis = [results_column,results_row];
        end
        
        %part b part 2

        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
        ThetaInDegrees = real(acosd(CosTheta));
        results_ang(results_column,results_row) = ThetaInDegrees;
        
        if final_result_ang > ThetaInDegrees
            final_result_ang = ThetaInDegrees;
            final_result_indicies_ang = [results_column,results_row];
        end
        
    end
end


results_dis
final_result_dis
final_result_indicies_dis

results_ang
final_result_ang
final_result_indicies_ang



document_frequency = zeros(1651,1);

for ind = 1:1651
    frequency = 0;
    for results_column = 1:10
        if W(results_column,ind)>0
            frequency = frequency + 1;
        end
        
    end
    document_frequency(ind) = sqrt(log(10/frequency));
end


for results_column = 2:10;
    for results_row = 1:(results_column-1);
        
        %part c part 1
        transposed_doc_freq = transpose(document_frequency);
        u_1 = W(results_column,:);
        u = dot((u_1/norm(u_1)),transposed_doc_freq);
        v_1 = W(results_row,:);
        v = dot((v_1/norm(v_1)),transposed_doc_freq);
        
        
        D = norm(u-v);
        results_dis(results_column,results_row) = D;
        
        if final_result_dis > D
            final_result_dis = D;
            final_result_indicies_dis = [results_column,results_row];
        end

        
    end
end

results_dis
final_result_dis
final_result_indicies_dis



