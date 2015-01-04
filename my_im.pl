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
  % TODO
  process_cost(Task, Core, Time), !.

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
  % TODO: trick the cuts
  findall(process_cost(Task, C, Time),
          process_cost(Task, C, Time), ProcCosts),
  minETList(ScheduleList, ProcCosts, process_cost(Task, Core, MinTime)),
  % important, only one solution here.
  process_cost(Task, Core, MinTime), !.

% required predicate to find the exact solution
find_optimal(result(Sol, ET)):-
  %findall(process_cost(A,B,C), process_cost(A,B,C), PCs),
  sort_tasks(PCs),
  extract_list(PCs, PCList),
  length(PCList, Len),
  print(Len),
  find_one(PCList, [], result(Sol, ET)).

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
minETList(ScheduleList, [process_cost(T1, C1, ET1),
                         process_cost(T2, C2, ET2)|Tail], Res):-

  add_to_schedule_list(process_cost(T1, C1, ET1),
                       ScheduleList,
                       NewScheduleList1),
  execution_time(solution(NewScheduleList1), TotalET1),
  add_to_schedule_list(process_cost(T2, C2, ET2),
                       ScheduleList,
                       NewScheduleList2),
  execution_time(solution(NewScheduleList2), TotalET2),


  ( (TotalET1 < TotalET2, find_dep(NewScheduleList1, T1)) ->
      minETList(ScheduleList, [process_cost(T1, C1, ET1)|Tail], Res)

  ;
    (find_dep(NewScheduleList2, T2) ->
        minETList(ScheduleList, [process_cost(T2, C2, ET2)|Tail], Res)
      ;
        minETList(ScheduleList, [process_cost(T1, C1, ET1)|Tail], Res)
    )

  ).

% gets the maximum processing cost in a list
minProcList([X], X) :- !.
minProcList([process_cost(T1, C1, Ti1),
             process_cost(T2, C2, Ti2)|Tail],
            Res):-
    ( Ti1 > Ti2 ->
        minProcList([process_cost(T2, C2, Ti2)|Tail], Res)
    ;
        minProcList([process_cost(T1, C1, Ti1)|Tail], Res)
    ).

% preprocessing: sorts task according to their execution time
sort_tasks(Res):-
  findall(T, task(T), Tasks),
  maplist(
    find_all(_),
  Tasks, MinPCs),
  predsort(compareDep, MinPCs, TmpList),
  predsort(compareTime, TmpList, Res).

%% helper function for maplist above
find_all(_, T, Res):-
  findall(process_cost(T, C, Ti), process_cost(T, C, Ti), Tasks),
  minProcList(Tasks, Res).

% compare function for predsort above

compareDep(<, process_cost(T1, _, _), process_cost(T2, _, _)):-
  depends_on(T1, T2, _), !.

compareDep(>, process_cost(T1, _, _), process_cost(T2, _, _)):-
  depends_on(T2, T1, _), !.

compareDep(<, process_cost(T1, _, _), process_cost(_, _, _)):-
  depends_on(T1, _, _), !.

compareDep(>, process_cost(_, _, _), process_cost(T2, _, _)):-
  depends_on(_, T2, _), !.

compareDep(>, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1<ET2, !.

compareDep(<, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1>=ET2, !.


% compare by time with respect to dependencies
compareTime(>, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  not(depends_on(T1, T2, _)),
  not(depends_on(T2, T1, _)),
  ET1<ET2, !.

compareTime(<, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  not(depends_on(T1, T2, _)),
  not(depends_on(T2, T1, _)),
  ET1>=ET2, !.

compareTime(<, process_cost(T1, _, _), process_cost(T2, _, _)):-
  depends_on(T1, T2, _).

compareTime(>, process_cost(T1, _, _), process_cost(T2, _, _)):-
  depends_on(T2, T1, _).

extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).
extract_list([], Tmp, Tmp).
extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).

% check dependencies
find_dep(_, DepTask):-
  not(depends_on(_, DepTask, _)), !.

find_dep(Schedules, DepTask):-
  depends_on(IndepTask, DepTask, _),
  get_schedule(Schedules, IndepTask, IndepSchdl),
  get_schedule(Schedules, DepTask, DepSchdl),
  tmp_exec_time(IndepSchdl, IndepTask, IndepET),
  tmp_exec_time(DepSchdl, DepTask, DepET),
  DepET > IndepET.


get_schedule([schedule(C, Ts)|_], Task, schedule(C, Ts)):-
  member(Task, Ts), !.

get_schedule([_|Schedules], Task, Res):-
  get_schedule(Schedules, Task, Res).

tmp_exec_time(Schedule, Task, Res):-
  tmp_exec_time(Schedule, Task, 0, Res).

tmp_exec_time(schedule(_, [T|_]), T, Tmp, Tmp):- !.

tmp_exec_time(schedule(C, [T|Ts]), Task, Tmp, Res):-
  process_cost(T, C, ET),
  Tmp1 is Tmp + ET,
  tmp_exec_time(schedule(C, Ts), Task, Tmp1, Res).
