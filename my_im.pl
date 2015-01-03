%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_one
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% helper for delete in lists
del(X,[X|Tail],Tail).
del(X,[Y|Tail],[Y|Tail1]):- del(X,Tail,Tail1).

% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).

list_min([L|Ls], Min) :-
    list_min(Ls, L, Min).

list_min([], Min, Min).
list_min([L|Ls], Min0, Min) :-
    Min1 is min(L, Min0),
    list_min(Ls, Min1, Min).


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

%  TODO: is Solution?
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
  schedule_execution_time([Schdl|Schdls], [], TotalET).

schedule_execution_time([], TotalETs, MaxTotalET):-
  maxList(TotalETs, MaxTotalET).

schedule_execution_time([schedule(Core,Ts)|Schdls], ETs, Res):-
  tasks_execution_time(Ts, Core, ET),
  append([ET], ETs, NewETs),
  schedule_execution_time(Schdls, NewETs, Res).

% Calculate execution time for schedule
tasks_execution_time([T|Ts], Core, TotalET):-
  process_cost(T, Core, TC),
  tasks_execution_time(Ts, Core, ET),
  TotalET is ET+TC.
tasks_execution_time([], _, TotalET):- TotalET is 0.

% get random processing cost
get_random_pc(Cores, Task, process_cost(Task, Core, Time)):-
  get_rndm_elem(Cores, Core),
  get_pc(process_cost(Task, Core, Time)), !.

find_heuristically(AllSolutions):-
  find_all_solutions(1000, [], AllSolutions).

find_all_solutions(0, SolList, SolList).

% create valid solution for the data above
find_all_solutions(Limit, SolList, ResSolList):-
  findall(T, task(T), AllTasks),
  find_one_solution(AllTasks, [], Sol),
  append(SolList, [Sol], NewSolList),
  Limit1 is Limit - 1,
  find_all_solutions(Limit1, NewSolList, ResSolList).


% create randomized solution
find_one_solution([Task|Tasks], ScheduleList, Solution):-
  findall(C, core(C), Cores),
  get_random_pc(Cores, Task, RndmPc),
  add_to_schedule_list(RndmPc, ScheduleList, NewScheduleList),
  find_one_solution(Tasks, NewScheduleList, Solution).

find_one_solution([], ScheduleList, ScheduleList).


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


% get list with random element (not in TakenElems)
get_rndm_elem(List, RndmElem):-
  length(List, Len),
  Len1 is Len - 1,
  random_between(0, Len1, RndmIndex),
  nth0(RndmIndex, List, RndmElem), !.

get_pc(ScheduleList, process_cost(Task, Core, MinTime)):-
  findall(process_cost(Task, C, Time),
          process_cost(Task, C, Time), ProcCosts),
  minETList(ScheduleList, ProcCosts, process_cost(Task, Core, MinTime)),
  % important, only one solution here.
  process_cost(Task, Core, MinTime), !.

% required predicate to find the exact solution
find_optimal(Res):-
  sort_tasks(PCs),
  extract_list(PCs, Tasks),
  find_one(Tasks, [], Res).

find_one([], ScheduleList, result(ScheduleList, ET)):-
  execution_time(solution(ScheduleList), ET), !.

find_one([Task|Tasks], ScheduleList, result(ResSchedule, ET)):-
  % binds only, if stated more solutions possible (?)
  get_pc(ScheduleList, process_cost(Task, Core, Time)),
  add_to_schedule_list(process_cost(Task, Core, Time),
                       ScheduleList,
                       NewScheduleList),
  find_one(Tasks, NewScheduleList, result(ResSchedule, ET)).

% checks candidates for a solution. chooses the task and cpu which
% increases the total execution time the least
minETList(_, [X], X) :- !.
minETList(ScheduleList, [ProcCost1, ProcCost2|Tail], Res):-
  add_to_schedule_list(ProcCost1,
                       ScheduleList,
                       NewScheduleList1),
  execution_time(solution(NewScheduleList1), ET1),
  add_to_schedule_list(ProcCost2,
                       ScheduleList,
                       NewScheduleList2),
  execution_time(solution(NewScheduleList2), ET2),
  ( ET1 > ET2 ->
      minETList(ScheduleList, [ProcCost2|Tail], Res)

  ;
      minETList(ScheduleList, [ProcCost1|Tail], Res)
  ).

% gets the maximum processing cost in a list
maxProcList([X], X) :- !.
maxProcList([process_cost(T1, C1, Ti1),
             process_cost(T2, C2, Ti2)|Tail],
            Res):-
    ( Ti1 > Ti2 ->
        maxProcList([process_cost(T1, C1, Ti1)|Tail], Res)
    ;
        maxProcList([process_cost(T2, C2, Ti2)|Tail], Res)
    ).

% sorts task according to their execution time
sort_tasks(Res):-
  findall(T, task(T), Tasks),
  maplist(
    find_all(_),
    Tasks, MaxTasks),
  predsort(compareTime, MaxTasks, Res).

% helper function for maplist above
find_all(_, T, Res):-
  findall(process_cost(T, C, Ti), process_cost(T, C, Ti), Tasks),
  maxProcList(Tasks, Res).

% compare function for predsort above
compareTime(Delta, process_cost(_, _, T1), process_cost(_, _, T2)):-
        compare(Delta, T2, T1).

extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).
extract_list([], Tmp, Tmp).
extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).

