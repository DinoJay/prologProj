%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_optimal
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


%get_random_pc(Cores, Task, process_cost(Task, Core, Time)):-
  %get_pc(process_cost(Task, Core, Time)).

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

find_optimal_pc(process_cost(Task, Core, Time),
                process_cost(Task, Core, Time)).

% workaround to bypass the cuts in 'batch_small_hetero'
get_pc(process_cost(Task, Core, Time)):- Core = c1,
  process_cost(Task, Core, Time).
get_pc(process_cost(Task, Core, Time)):- Core = c2,
  process_cost(Task, Core, Time).
get_pc(process_cost(Task, Core, Time)):- Core = c3,
  process_cost(Task, Core, Time).
get_pc(process_cost(Task, Core, Time)):- Core = c4,
  process_cost(Task, Core, Time).
get_pc(process_cost(Task, Core, Time)):-
  process_cost(Task, Core, Time).

% Testing
find_one(result(Res, ET)):-
  findall(T, task(T), Tasks),
  find_optimal(Tasks, [], result(Res, ET)).

% adds task schedule according the process_cost.
% Important, here the execution tree gets wider
find_optimal_pc(Task, ScheduleList, NewSol):-
  % crucial point, more solutions are possible for the following line
  %process_cost(Task, Core, Time),
  process_cost(Task, Core, Time),
  add_to_schedule_list(process_cost(Task, Core, Time), ScheduleList,
  NewSol).

find_optimal([], Schedules, result(Schedules, ET)):-
  execution_time(solution(Schedules), ET).

find_optimal([Task|Tasks], ScheduleList, result(ResSchedule, ET)):-
  % binds only, if stated more solutions possible (?)
  find_optimal_pc(Task, ScheduleList, NewScheduleList),
  find_optimal(Tasks, NewScheduleList, result(ResSchedule, ET)).

find_optimal(Res):-
  findall(T, task(T), Tasks),
  findall(Sol, find_optimal(Tasks, [], Sol), List),
  minList(List, Res).
  %find_min_ec(SolList, Res),
  %execution_time(solution(Res), ET).


minList([X], X) :- !.

minList([result(Schedule, ET),
                   result(Schedule1, ET1)|Tail],
                   Res):-
    ( ET > ET1 ->
        minList([result(Schedule1, ET1)|Tail], Res)
    ;
        minList([result(Schedule, ET)|Tail], Res)
    ).

find_min_ec([Sol|Sols], CurMinSol, MinSol):-
  execution_time(solution(CurMinSol), CurMin),
  execution_time(solution(Sol), ET),
  ( CurMin > ET ->
      CurMinSol1 = Sol;
      CurMinSol1 = CurMinSol
  ),
  find_min_ec(Sols, CurMinSol1, MinSol).

find_min_ec([Sol|Sols], Min):-
  find_min_ec(Sols, Sol, Min).

% find min process_cost() in list
