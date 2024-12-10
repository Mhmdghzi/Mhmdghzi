clear all
clc
NC = [200];             % number of circles
r = 0.235;           % radius of circle
%% For BS1
th = linspace(0, 2*pi);            % Angles for circle circumference
x_circle = r * cos(th);            % x-coordinates for circle circumference
y_circle = r * sin(th);            % y-coordinates for circle circumference
steps = [0:1:170];     % Steps for x-coordinate of the second point
run_times=5;
 blockers_mean=zeros(run_times,length(steps));       %blockers matrix
  blockers_mean2=zeros(run_times,length(steps));


  % mkdir('C:\Users\slo0o\Downloads\LinePath3\');
 %save_dir = 'C:\Users\slo0o\OneDrive\سطح المكتب\project reults\LinePath results\markov fix attempt\New folder\';

start_time = datetime('now');   % just to see how long it will take the code to run

width = 170; % Width of the rectangle
height = 15; % Height of the rectangle


% for NB=1:length(NC)
N=NC;
  filename =['linepath',num2str(N)];
for  b=1:run_times
  
C = [width * rand(N, 1), height * rand(N, 1)];

    disp(b)
%     start_point = 0; % Choose the start point
% end_point = 40;   % Choose the end point
% C = start_point + (end_point - start_point) * rand(N, 2);
 blockers = zeros(1, length(steps));            % Initialize the counter for blockers
for step = 1:numel(steps)       
     
    for i = 1:N
        xc = C(i, 1) + x_circle;    % x-coordinates for circle
        yc = C(i, 2) + y_circle;    % y-coordinates for circle
    
        % Define the line segment with varying x-coordinate
        x = [steps(step) 85];   % Increase x-coordinate by 0.2 in each iteration
        y = [0 15];                  % Fixed y-coordinate
        theta=atan(x(1)/y(2));              % calculating the length of the BS
        r=y(2)/theta;     
        
       % axis equal
       %hold on
       %plot(xc, yc);
        % Find intersection points between circle and line segment
        [xi, yi] = polyxpoly(x, y, xc, yc);
          %Plot intersection points and line segment
        %  mapshow(xi, yi, 'DisplayType', 'point', 'Marker', 'o');
        %  mapshow(x, y, 'Marker', '+');
        % shows the plots
        if sqrt((x(2)-C(i, 1))^2+(y(2)-C(i, 2))^2)>=r     % ensuring blockers in the no-blockers area dont effect
        if ~isempty(xi)
            blockers(step) = blockers(step) + 1;           % number of blockers
             blockers_mean(b,step)=blockers(step);       %the matrix of the blocker                         
        end
        end
    end
    
   
end


end






%% 6markov chain


number_of_sections=12;



%  this loop will divied the blockers_mean matrix into six matrixes
 blockers_mean_divisions=cell(1,number_of_sections);
 division_factor=floor(length(steps)/number_of_sections);
 starts=1;
addition=division_factor;

transition_matrix=cell(1,number_of_sections);     % making transition matrix to each division
transition_probability_matrix=cell(1,number_of_sections);
markov_division=cell(1,number_of_sections);
markov=zeros(length(blockers_mean(:,1)),length(blockers_mean(1,:)));
for x=1:number_of_sections


blockers_mean_divisions{x}=blockers_mean(:,starts:addition);




addition=division_factor+addition;

starts=starts+division_factor;










%% making transition matrix



  states=length(unique(blockers_mean_divisions{x}));

        states_occurnces=zeros(1,states);

        for i=0:(length(steps)-1)    %making the loop run as many times as the steps to cover all columns.

        states_occurnces(i+1)=sum(sum(blockers_mean_divisions{x}==i));  %counting the  state occurnces for all the matrix

        end



 
  % making the matrix for trasition from a number to a speecific number  {n(i-j)}

 transition_matrix{x}=zeros(states);
  for b=0:states-1   % this loop for the start value for current state, and look on that value to cover the entire matrix
 for  a=0:run_times-1 % this loop for cover all the raws in blockers_mean
 for i=1:length(blockers_mean_divisions{x}(1,:))-1     % this loop to cover all the columns in blockers_mean
     current_state=blockers_mean_divisions{x}(a+1,i);
     next_state=blockers_mean_divisions{x}(a+1,i+1);
     for q=0:states-1          % this loop to check the next state with all possibel values before going to the next psition
     if  current_state==b && next_state==q    
         transition_matrix{x}(b+1,q+1)=transition_matrix{x}(b+1,q+1)+1;  

     end
     end

 end
 end
 
  end
















  % this loop to check for the last value in each row and replicate it, 
 for a=0:states-1     % to check the vlaues from 0 to 5 if they are in the last column from each row, since these are numbers that occur
     for s=1:run_times   % checking in all the rows in matrix blockers_mean
             if blockers_mean_divisions{x}(s,end)==a
           states_occurnces(a+1)=states_occurnces(a+1)+1; 

             transition_matrix{x}(a+1,a+1)=transition_matrix{x}(a+1,a+1)+1;   

             end
     end
 end

 % this is the final transition probability matrix, divding each element
 % from sum_transition_matrix with a whole row of transition_matrix
 
 transition_probability_matrix{x}=zeros(states,states);
 sum_transition_matrix=sum(transition_matrix{x},2)';

 %this loop to check if there is a line of zeros in (transition_matrix and)
 %remove it, so there won't be NaN case in the (transition_probability_matrix)
 % for a=1:length(transition_matrix{x}(:,1))
 % 
 %     if transition_matrix{x}(a,:)==0
 %        transition_matrix{x}(a,:)=[];
 % 
 %     end
 % end
 % 
TR=transition_probability_matrix{x};

for a=1:length(transition_matrix{x}(:,1))
 for i=1:length(transition_matrix{x}(1,:))

  transition_probability_matrix{x}(a,i) =transition_matrix{x}(a,i)/sum_transition_matrix(a);

 end
end



%%  end of transition matrix
%save('C:\Users\slo0o\Downloads\LinePath\transition_probability_linePath.mat','transition_probability_matrix');   %this line to save the variable
%"transition_probability_matrix" 








%%  markov model



markov_division{x}=zeros(length(blockers_mean_divisions{x}(:,1)),length(blockers_mean_divisions{x}(1,:)));     % generate the markov matrix that is the exact size as the blockers_mean mtrix
temp_transition_probability_matrix=transition_probability_matrix;


for a=1:length(blockers_mean_divisions{x}(:,1))    % running throught all the rows blockers_mean
%initial_state=randsample(unique(blockers_mean_divisions{x}),1);     %statrting the initial state of every row 
initial_state=round(mean(blockers_mean_divisions{x}(a,:)));     %statrting the initial state of every row 
% if x>1
%     initial_state=blockers_mean_divisions{x-1}(a,1);
% end

markov_division{x}(a,1)=initial_state;    % starting the initial_state of of each run for the markov
for i=1:length(blockers_mean_divisions{x}(1,:))-1   % running through all the columns blockers_mean
next_state=rand(1);      %  the prapobility value
    v=initial_state;
% if next_state<=transition_probability_matrix(initial_state+1,m)
m=1;     %a value to statrt in the first value of each row in transition_probability_matrix
if initial_state==length(transition_probability_matrix{x}(:,1))
   initial_state= initial_state-1;
end
current_state=transition_probability_matrix{x}(initial_state+1,m);

while next_state>=current_state

current_state=current_state+transition_probability_matrix{x}(initial_state+1,m+1);
m=m+1;             % "m" here have the value of the next position for the matrix markov
% disp('next_state ')
% disp(next_state);
% disp('current_state ')
% disp(current_state);
end



     markov_division{x}(a,i+1)=m-1;
end

end


end

DV=division_factor;
CV=1;
for i=1:number_of_sections

    markov(:,CV:DV)=  markov_division{i};

CV=DV+1;
DV=DV+division_factor;

end


mean_absolute_error = mean(abs(mean(markov) - mean(blockers_mean)));

disp(['mean_absolute_error=:',num2str(mean_absolute_error)])


markov(:,1)=[];
blockers_mean(:,1)=[];
a=mean(markov);
b=mean(blockers_mean);
figure;                % Create a new figure window
plot(a, 'b-');         % Plot 'a' with a solid blue line
hold on;               % Retain the current plot so we can add another plot

plot(b, 'r--');        % Plot 'b' with a dashed red line

% Create a dynamic title using the variables for steps, runs, and blockers
the_title = ['strsight Line Path, Steps: ', num2str(length(steps)), ', Runs: ', num2str(run_times), ', Blockers: ', num2str(N), ',sections:',num2str(number_of_sections) ];
title(the_title);      % Add title to the plot

xlabel('Steps');       % Label for the x-axis
ylabel('Average Blockers');   % Label for the y-axis

hold off; 

% figure;
% plot(b,'r--')
% the_title = ['markov prediction line Path, Steps: ', num2str(length(steps)), ', Runs: ', num2str(run_times),',blockers:',num2str(N)];
% title(the_title); % Add title
% xlabel('steps'); % Add x-axis label
% ylabel('average blocckers'); % Add y-axis label
% 






% Construct the filename using the variables in your code
%filename = ['blockers_mean_runs:', num2str(run_times), '_steps', 5num2str(length(steps)), '_BL', num2str(NC), '.mat'];

% Save the 'blockers_mean' matrix with the constructed filename
%save(fullfile('C:\Users\slo0o\Downloads\matlab\line path\', filename), 'blockers_mean');



























% markov(:,1)=[];
% blockers_mean(:,1)=[];
% a=mean(markov);
% b=mean(blockers_mean);
% figure; 
% plot(a,'b-')
% hold on
% plot(b,'r--')
% hold off; % Release current plot
% the_title = ['line Path, Steps: ', num2str(length(steps)), ', Runs: ', num2str(run_times),',blockers:',num2str(N)];
% title(the_title); % Add title
% xlabel('steps'); % Add x-axis label
% ylabel('average blocckers'); % Add y-axis label
% legend('markov', 'actual data ', 'Location', 'northwest'); % Add legend in the top left corner


 % saveas(gcf, fullfile(save_dir, filename), 'jpg'); % fullfile to combine path and filename
    % close(gcf);







end_time = datetime('now');
%Calculate the elapsed time
elapsed_time = end_time - start_time;
disp(['Elapsed time: ', char(elapsed_time)]);



%% end of markov






% 
% for i=0:length(blockers_mean(1,:))-1
% 
% 
% if blockers_mean(1,1)==i
% markov(i,1)=transition_probability_matrix(i+1,i);
% 
% random=rand(1,1);
% next_state=transition_probability_matrix(i+1,:);
% 
% for b=1:length(transition_probability_matrix(1,:))-1
% 
% if random <= transition_probability_matrix(1,b)
% 
% %next_state=
% 
% 
% 
% end
% end
% end
% end        
