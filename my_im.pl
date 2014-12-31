%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% create solution, find_optimal
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% helper for delete in lists
del(X,[X|Tail],Tail).
del(X,[Y|Tail],[Y|Tail1]):- del(X,Tail,Tail1).

% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).

% find min process_cost() in list
minList([X], X) :- !.
minList([process_cost(T,C,PC),
                   process_cost(T1,C1, PC1)|Tail],
                   Res):-
    ( PC > PC1 ->
        minList([process_cost(T1,C1,PC1)|Tail], Res)
    ;
        minList([process_cost(T,C, PC)|Tail], Res)
    ).


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

schedule_execution_time([], TotalETs, MaxTotalET):-
  maxList(TotalETs, MaxTotalET).

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
create_solution(Solution):-
  create_solution([], [], [], Solution).

create_solution(ScheduleList, _, _, ScheduleList):-
  isSolution(solution(ScheduleList)), !.

% create randomized solution
create_solution(ScheduleList, TakenCores, TakenTasks, Solution):-

  create_schedule(TakenCores, TakenTasks, NewTakenCores, NewTakenTasks,
                  schedule(Core, TaskSet)),

  append([schedule(Core, TaskSet)], ScheduleList, NewScheduleList),
  create_solution(NewScheduleList, NewTakenCores, NewTakenTasks, Solution).


% create randomized schedule
create_schedule(TakenCores, TakenTasks, NewTakenCores, NewTakenTasks,
                schedule(Core, TaskSet)):-

  findall(C, core(C), Cores),
  get_rndm_elem(Cores, TakenCores, Core),

  findall(T, task(T), AllTasks),
  gen_rndm_list(AllTasks, Core, TakenTasks, TaskSet),

  append(TakenCores, [Core], NewTakenCores),
  append(TakenTasks, TaskSet, NewTakenTasks).


% get list with random elements (not in TakenElems)
gen_rndm_list(List, Core, TakenElems, RES):-
  length(List, Len),
  Len1 is Len - 1,
  random_between(1, Len1, RndmLen),
  gen_rndm_list(RndmLen, List, Core, TakenElems, RES), !.

gen_rndm_list(0, _, _, []).
gen_rndm_list(_, L1, L2, []):- subtract(L1,L2, []), !.

gen_rndm_list(Len, Core, List, TakenElems, RES):-
  get_rndm_elem(List, Core, TakenElems, RndmTask),
  Len1 is Len - 1,
  append(TakenElems, [RndmTask], NewTakenElems),
  gen_rndm_list(Len1, Core, List, NewTakenElems, SubRes),
  append([RndmTask], SubRes, RES),
  !.

% for testing purposes
%add_to_schedules(NewSol):-
  %create_solution(Sol),
  %find_optimal_pc(t2, PC),
  %% add new task to schedule in possible solution
  %add_to_schedule_list(PC, Sol, NewSol).

add_to_schedule_list(process_cost(Task, Core, Cost), Sol, NewSol):-
  maplist(add(process_cost(Task, Core, Cost)), Sol, NewSolTmp),
  (Sol = NewSolTmp ->
    append(NewSolTmp, [schedule(Core, [Task])], NewSol);
    append(NewSolTmp, [], NewSol)
  ).

add(process_cost(Task, Core, _), schedule(Core, Tasks),
  schedule(Core, NewTasks)):-
      append(Tasks, [Task], NewTasks),
      !.

add(process_cost(_, _, _), schedule(Core1, Tasks),
  schedule(Core1, Tasks)).

copy(A, A).

gen_random(Core, Remainder, Res):-
  random_between(0, Len1, RndmIndex),
  length(Remainder, Len),
  Len1 is Len - 1,
  nth0(RndmIndex, Remainder, Task),
  (process_cost(Task, Core, _) ->
    copy(Task, Res);
    gen_random(Core, Remainder, Res)
  ).

% get list with random element (not in TakenElems)
get_rndm_elem(List, TakenElements, RndmElem):-
  subtract(List, TakenElements, Remainder),
  length(Remainder, Len),
  Len1 is Len - 1,
  random_between(0, Len1, RndmIndex),
  nth0(RndmIndex, Remainder, RndmElem), !.

get_rndm_elem(List, Core, TakenElements, RndmElem):-
  subtract(List, TakenElements, Remainder),
  gen_random(Core, Remainder, RndmElem).

find_optimal_pc(process_cost(Task, Core, Time),
                process_cost(Task, Core, Time)).

% workaround to bypass the cuts in 'batch_small_hetero'
get_pc(process_cost(_, Core, _)):- Core = c1.
get_pc(process_cost(_, Core, _)):- Core = c2.
get_pc(process_cost(_, Core, _)):- Core = c3.
get_pc(process_cost(_, Core, _)):- Core = c4.
get_pc(process_cost(_, _, _)).

% adds task schedule according the process_cost.
% Important, here the execution tree gets wider
find_optimal_pc(Task, ScheduleList, NewSol):-
  get_pc(process_cost(Task, Core, Cost)),
  add_to_schedule_list(process_cost(Task, Core, Cost), ScheduleList,
  NewSol).

find_optimal([], Schedules, Schedules).

find_optimal([Task|Tasks], ScheduleList, Res):-
  % binds only, if stated more solutions possible (?)
  find_optimal_pc(Task, ScheduleList, NewScheduleList),
  find_optimal(Tasks, NewScheduleList, Res).
  %isSolution(isSolution(C)).
  % TODO: get the schedule with the minimum execution time

find_optimal(Res):-
  findall(T, task(T), Tasks),
  find_optimal(Tasks, [], Res).


