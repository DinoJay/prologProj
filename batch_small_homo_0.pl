%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          batch_small_homo.pl          %
%                                       %
%    Scheduling of a small batch of     %
%   independent tasks on a homogeoneous %
%                system                 %
%                                       %
%       Declarative Programming         %
%              2014-2015                %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%%   Cores  %%%
%%%%%%%%%%%%%%%%
% The cores of the machine
% core(Id): a core with unique identifier 'Id'.
core(c1).
core(c2).
core(c3).
core(c4).

%%%%%%%%%%%%%%%%%%
%%%   Tasks    %%%
%%%%%%%%%%%%%%%%%%
% The tasks the application is made up of, which are to be scheduled.
% task(Id): a task with unique identifier 'Id'.
task(t1).
task(t2).
task(t3).
task(t4).
task(t5).
task(t6).
task(t7).

%%%%%%%%%%%%%%%%%%
%% Dependencies %%
%%%%%%%%%%%%%%%%%%
% The execution order dependencies between tasks
% depends_on(Ta,Tb,Data): before task 'Ta' can be executed,
% task 'Tb' must have been executed and thereafter 'Data' megabytes of data (result of/produced by 'Tb') must have been moved from the processor that executed 'Tb' to the processor that will execute 'Ta'.

%In this benchmark there are no dependencies between tasks.
depends_on(_,_,_) :- fail.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Processing Costs   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies how long the processing of each task takes on each core.
% process_cost(T,C,Time): It takes 'Time' ms to execute task 'T' on core 'C'.

%In this benchmark every core can execute a given task as fast (i.e. homogeneous system)
process_cost(t1,C,100) :- core(C).
process_cost(t2,C,20) :- core(C).
process_cost(t3,C,30) :- core(C).
process_cost(t4,C,40) :- core(C).
process_cost(t5,C,60) :- core(C).
process_cost(t6,C,70) :- core(C).
process_cost(t7,C,80) :- core(C).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Channel Properties  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies the properties of the communication channel between each of the cores
% channel(Ca,Cb,Latency,Bandwidth): The channel to communicate from core 'Ca' to core 'Cb' has a latency 'Latency' and bandwidth 'Bandwidth'.
% Note that sending 'X' megabytes of data, over a channel, takes Latency + X/Bandwidth ms.

%In this benchmark (without task dependencies), no inter core communication is required, the channel properties therefore do not matter.
channel(C1,C2,0,1) :- core(C1), core(C2).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  My Implementation %%%
% quite sketchy, but first importan function, create solution,
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% helper for delete in lists
del(X,[X|Tail],Tail).
del(X,[Y|Tail],[Y|Tail1]):- del(X,Tail,Tail1).

% find max in list
maxList([A],A).
maxList([A|List],Max):-
 maxList(List,Max1),
 (A>=Max1, Max=A; A<Max1, Max=Max1).

% valid schedule --> valid core and valid tasks
isSchedule(schedule(C, [T|Ts])):- core(C), isTaskSet(T, Ts).
isTaskSet(T, [To|Ts]):- task(T), isTaskSet(To, Ts).
isTaskSet(T, []):- task(T).

% mark tasks of a schedule in a given tasks set
markTasks(schedule(_, Tasks), AllTasks, ReturnTasks):-
  markTasks(Tasks, AllTasks, ReturnTasks).

markTasks([T|Ts], AllTasks, ReturnTasks):-
  del(T, AllTasks, NewAllTasks),
  markTasks(Ts, NewAllTasks, ReturnTasks).

markTasks([], NewAllTasks, NewAllTasks).

% is Solution?
isSolution(solution(Schedules)):-
  findall(T, task(T), AllTasks),
  isSolution(Schedules, AllTasks).

isSolution([Schdl|Schdls], AllTasks):-
  isSchedule(Schdl),
  markTasks(Schdl, AllTasks, NotProcTasks),
  isSolution(Schdls, NotProcTasks).

isSolution([], []).
isSolution([], _):- false.


% Calculate execution time for solution
execution_time(solution([Schdl|Schdls]), TotalET):-
  isSolution(solution([Schdl|Schdls])),
  schedule_execution_time([Schdl|Schdls], [], TotalET).

schedule_execution_time([], Res, MaxRes):- maxList(Res, MaxRes).

schedule_execution_time([schedule(_,Ts)|Schdls], ETs, Res):-
  tasks_execution_time(Ts, ET),
  append([ET], ETs, NewETs),
  schedule_execution_time(Schdls, NewETs, Res).

% Calculate execution time for schedule
tasks_execution_time([T|Ts], TotalET):-
  process_cost(T, _, TC),
  tasks_execution_time(Ts, ET),
  TotalET is ET+TC.

tasks_execution_time([], TotalET):- TotalET is 0.


% create valid solution for the data above
create_solution(TotalET):-
  create_solution([], [], [], Solution),
  schedule_execution_time(Solution, [], TotalET).


create_solution(ScheduleList, _, _, ScheduleList):-
  isSolution(solution(ScheduleList)), !.

create_solution(ScheduleList, TakenCores, TakenTasks, Solution):-
  create_schedule(TakenCores, TakenTasks, schedule(Core, TaskSet)),
  append(TakenCores, [Core], NewCores),
  append(TakenTasks, TaskSet, NewTasks),
  list_to_set(NewTasks, NewTakenTasks),
  list_to_set(NewCores, NewTakenCores),
  append([schedule(Core, TaskSet)], ScheduleList, NewScheduleList),
  create_solution(NewScheduleList, NewTakenCores, NewTakenTasks, Solution).


create_schedule(TakenCores, TakenTasks, schedule(RndmCore, TaskSet)):-
  findall(C, core(C), Cores),
  get_rndm_elem(Cores, TakenCores, RndmCore),

  findall(T, task(T), AllTasks),
  gen_rndm_list(AllTasks, TakenTasks, TaskSet).

% get list with random elements (not in TakenElems)
gen_rndm_list(List, TakenElems, RES):-
  length(List, Len),
  Len1 is Len - 1,
  random_between(1, Len1, RndmLen),
  gen_rndm_list(RndmLen, List, TakenElems, RES), !.

gen_rndm_list(0, _, _, []).
gen_rndm_list(_, L1, L2, []):- subtract(L1,L2, []), !.

gen_rndm_list(Len, List, TakenElems, RES):-
  get_rndm_elem(List, TakenElems, RndmTask),
  Len1 is Len - 1,
  append(TakenElems, [RndmTask], NewTakenElems),
  gen_rndm_list(Len1, List, NewTakenElems, SubRes),
  append([RndmTask], SubRes, RES),
  !.


% get list with random element (not in TakenElems)
get_rndm_elem(List, TakenElements, RndmElem):-
  subtract(List, TakenElements, Rem),
  length(Rem, Len),
  Len1 is Len - 1,
  random_between(0, Len1, RndmIndex),
  nth0(RndmIndex, Rem, RndmElem), !.

% TODO:find optimal solution
%find_optimal(C):-
  %isSolution(isSolution(C)).

