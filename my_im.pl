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
  create_schedule(TakenCores, TakenTasks, schedule(Core, TaskSet)),
  append(TakenCores, [Core], NewCores),
  append(TakenTasks, TaskSet, NewTasks),
  list_to_set(NewTasks, NewTakenTasks),
  list_to_set(NewCores, NewTakenCores),
  append([schedule(Core, TaskSet)], ScheduleList, NewScheduleList),
  create_solution(NewScheduleList, NewTakenCores, NewTakenTasks, Solution).

% create randomized schedule
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

add_to_schedules(NewSol):-
  create_solution(Sol),
  find_optimal_PC(t2, PC),
  % add new task to schedule in possible solution
  add_to_schedule_list(PC, Sol, NewSol).

add_to_schedule_list(process_cost(Task, Core, Cost), Sol, NewSol):-
  maplist(add(process_cost(Task, Core, Cost)), Sol, NewSolTmp),
  (Sol = NewSolTmp ->
    append(NewSolTmp, [schedule(Core, [Task])], NewSol);
    append(NewSolTmp, [], NewSol)
  ).


%add_to_schedules(PC, [], NewSol):-
  %maplist(add(PC), [], NewSol).

add(process_cost(Task, Core, _), schedule(Core, Tasks),
  schedule(Core, NewTasks)):-
      append(Tasks, [Task], NewTasks),
      !.

add(process_cost(_, _, _), schedule(Core1, Tasks),
  schedule(Core1, Tasks)).


% get list with random element (not in TakenElems)
get_rndm_elem(List, TakenElements, RndmElem):-
  subtract(List, TakenElements, Rem),
  length(Rem, Len),
  Len1 is Len - 1,
  random_between(0, Len1, RndmIndex),
  nth0(RndmIndex, Rem, RndmElem), !.


find_optimal_PC(Task, MinPC):-
  findall(process_cost(Task, Core, Cost),
          process_cost(Task, Core, Cost), PCTasks),
  minList(PCTasks, MinPC).

find_optimal([], Schedules, Schedules).

find_optimal([Task|Tasks], Schedules, Res):-
  find_optimal_PC(Task, MinPC),
  add_to_schedule_list(MinPC, Schedules, NewSchedules),
  find_optimal(Tasks, NewSchedules, Res).
  %isSolution(isSolution(C)).

find_optimal(Res):-
  findall(T, task(T), Tasks),
  find_optimal(Tasks, [], Res).



