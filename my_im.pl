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

minList([L|Ls], Min) :-
    minList(Ls, L, Min).

minList([], Min, Min).
minList([L|Ls], Min0, Min) :-
    Min1 is min(L, Min0),
    minList(Ls, Min1, Min).


% valid schedule --> valid core and valid tasks
isSchedule(schedule(C, [T|Ts])):- core(C), isTaskSet(T, Ts).
isTaskSet(T, [To|Ts]):- task(T), isTaskSet(To, Ts).
isTaskSet(T, []):- task(T).


%% mark tasks of a schedule in a given tasks set
%markTasks(schedule(_, Tasks), AllTasks, ReturnTasks):-
  %markTasks(Tasks, AllTasks, ReturnTasks).

%markTasks([T|Ts], AllTasks, ReturnTasks):-
  %del(T, AllTasks, NewAllTasks),
  %markTasks(Ts, NewAllTasks, ReturnTasks).

%markTasks([], NewAllTasks, NewAllTasks).

%  TODO: is Solution?
isSolution(solution(Schedules)):-
  findall(T, task(T), AllTasks),
  length(AllTasks, Len),
  isSolution(Schedules, 0, Num),
  Len = Num.

isSolution([schedule(C, Ts)|Schdls], Counter, Res):-
  length(Ts, NumTs),
  Counter1 is Counter + NumTs,
  isSolution(Schdls, Counter1, Res).

isSolution([], Counter, Counter).


% get random processing cost
get_random_pc(Cores, Task, process_cost(Task, Core, Time)):-
  get_rndm_elem(Cores, Core),
  process_cost(Task, Core, Time), !.

find_heuristically(Limit, AllSolutions):-
  find_heuristically(Limit, [], AllSolutions).

find_heuristically(0, SolList, Res):- minList(SolList, Res).

% find solutions corresponding to the given limit
find_heuristically(Limit, SolList, ResSolList):-
  findall(T, task(T), AllTasks),
  find_one_solution(AllTasks, [], Sol),
  execution_time(Sol, ET),
  append(SolList, [ET], NewSolList),
  Limit1 is Limit - 1,
  find_heuristically(Limit1, NewSolList, ResSolList), !.


% create randomized solution
find_one_solution([Task|Tasks], ScheduleList, Solution):-
  findall(C, core(C), Cores),
  get_random_pc(Cores, Task, RndmPc),
  add_to_schedule_list(RndmPc, ScheduleList, NewScheduleList),
  find_one_solution(Tasks, NewScheduleList, Solution).

find_one_solution([], ScheduleList, ScheduleList).

% add process_cost(T, C, ET) to the corresponding schedule
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
  depSort(TaskList),
  find_one(TaskList, [], result(Sol, ET)), !.

find_one([], ScheduleList, result(ScheduleList, ET)):-
  execution_time(ScheduleList, ET), !.

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
  execution_time(NewScheduleList1, TotalET1),
  add_to_schedule_list(process_cost(T2, C2, ET2),
                       ScheduleList,
                       NewScheduleList2),
  execution_time(NewScheduleList2, TotalET2),


  ( (TotalET1 < TotalET2) ->
      minETList(ScheduleList, [process_cost(T1, C1, ET1)|Tail], Res)

  ;
      minETList(ScheduleList, [process_cost(T2, C2, ET2)|Tail], Res)
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


maxETList([X], X) :- !.
maxETList([process_cost(T1, C1, Ti1),
             process_cost(T2, C2, Ti2)|Tail],
            Res):-
    ( Ti1 > Ti2 ->
        minProcList([process_cost(T2, C2, Ti2)|Tail], Res)
    ;
        minProcList([process_cost(T1, C1, Ti1)|Tail], Res)
    ).

% preprocessing: sorts task according to their execution time
sort_tasks(Tasks, Res):-
  maplist(
    find_all(_),
  Tasks, PCs),
  flatten(PCs, FlattedPCs),
  predsort(compareTime, FlattedPCs, Res).

%% helper function for maplist above
find_all(_, T, Res):-
  findall(process_cost(T, C, Ti), process_cost(T, C, Ti), Res).
  minProcList(Tasks, Res).

getPCs(Paths, Res):- getPCs(Paths, [], Res).
getPCs([], Tmp, Res):- list_to_set(Tmp, Res).
getPCs([path(Ts)|Tail], Tmp, Res):-
  append(Tmp, Ts, Tmp1),
  getPCs(Tail, Tmp1, Res).

% sorting based on dependencies, i.e. topological sorting
depSort(Sorted):-
  findall(T, make_edge(T), Ts),
  (Ts = [] ->
    fail
  ;
    flatten(Ts, Fl),
    vertices_edges_to_ugraph([], Fl, L),
    top_sort(L, Sorted), !
  ).

depSort(Res):-
  findall(T, task(T), Ts),
  sort_tasks(Ts, SortedPCs),
  extract_list(SortedPCs, SortedTs),
  list_to_set(SortedTs, Res).


make_edge([EnablingTask-DepTask]):- depends_on(DepTask, EnablingTask, _).

% compare by time with respect to dependencies
compareTime(>, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1<ET2.

compareTime(<, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1>=ET2.

extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).
extract_list([], Tmp, Tmp).
extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).

is_sink(AllTasks, Task):-
  findall(T, depends_on(T, Task, _), DepTasks),
  subtract(AllTasks, DepTasks, DiffTasks),
  DiffTasks = AllTasks.


find_sink(AllTasks, Res):-
  find_sink(AllTasks, AllTasks, Res).

find_sink(AllTasks, [T|_], T):-
  is_sink(AllTasks, T).

find_sink(AllTasks, [_|Ts], Res):-
  find_sink(AllTasks, Ts, Res).

execution_time(Schedules, Res):-
  execution_time(Schedules, Schedules, ETs),
  flatten(ETs, FlattedETs),
  maxList(FlattedETs, Res),
  !.

execution_time(ScheduleList, ScheduleList, Res):-
  execution_time(ScheduleList, ScheduleList, [], Res).

execution_time(ScheduleList, [], Tmp, Tmp).

execution_time(ScheduleList, [schedule(C, Ts)|Tail], Tmp, Res):-
  % get all paths to the root
  findall(ET,
    execution_time(ScheduleList, schedule(C, Ts), ET),
  NestedPaths),
  flatten(NestedPaths, Paths),
  % get longest (most expensive) path to the root
  maxPCs(Paths, tuple(MaxPath, MaxET)),
  % get execution times for schedule
  getPCWrapper([schedule(C, Ts)], Ts, SchedulePCs),
  %getPCs(Paths, PathPCs),
  subtract(SchedulePCs, MaxPath, DiffPCs),
  sumPCs(DiffPCs, DiffPCsSum),
  Max is DiffPCsSum + MaxET,
  append(Tmp, [Max], Tmp1),
  execution_time(ScheduleList, Tail, Tmp1, Res), !.

execution_time(ScheduleList, schedule(C, Ts), Res):-
  % Several Solutions
  find_sink(Ts, Sink),
  findall(path(Path),
    search_df(ScheduleList, schedule(C, Ts), [[Sink]], Path),
  Res).


% depth first search
search_df(ScheduleList, Schedule, [Current|_], Res):-
  children(Current,Children),
  Children = [],
  getPCWrapper(ScheduleList, Current, Res).

search_df(ScheduleList, Schedule, [Current|Rest], Res):-
  children(Current, Children),
  append(Children, Rest, NewAgenda),
  search_df(ScheduleList, Schedule, NewAgenda, Res).

children([Node|RestOfPath], Children):-
 findall([Child, Node|RestOfPath], depends_on(Node, Child, _), Children),
 !.


% get process_costs for schedule
getPCWrapper(Schedules, Tasks, Res):-
  getPCWrapper(Schedules, Tasks, [], Res).

getPCWrapper(Schedules, [], Tmp, Tmp).

getPCWrapper(Schedules, [T|Ts], Tmp, Res):-
  getPC(Schedules, T, PC),
  append([PC], Tmp, Tmp1),
  getPCWrapper(Schedules, Ts, Tmp1, Res).

getPC([Schedule|Tail], Task, Res):-
  getPC(Schedule, Task, Res).

getPC([_|Tail], Task, Res):-
  getPC(Tail, Task, Res).

getPC(schedule(C, Ts), Task, process_cost(Task, C, ET)):-
  member(Task, Ts),
  process_cost(Task, C, ET), !.

getPC(schedule(C, [_|Ts]), Task, Res):-
  getPC(schedule(C, Ts), Task, Res).

getPC(schedule(_, []), _, _):- fail.


maxPCs([path(PCs)|Tail], Res):-
    sumPCs(PCs, Max),
    maxPCs(Tail, tuple(PCs, Max), Res).

maxPCs([], MaxTuple, MaxTuple).
maxPCs([path(PCs)|Tail], tuple(PCs0, Max0), Res):-
  sumPCs(PCs, Sum),
  (Sum > Max0 ->
    maxPCs(Tail, tuple(PCs, Sum), Res)
  ;
    maxPCs(Tail, tuple(PCs0, Max0), Res)
  ).


sumPCs(PCs, Res):-
  sumPCs(PCs, 0, Res).

sumPCs([], Tmp, Tmp).

sumPCs([process_cost(T0, C0, ET0),
        process_cost(T1, C1, ET1)|PCs], Tmp, Res):-
  % Note that sending 'X' megabytes of data, over a channel,
  % takes Latency + X/Bandwidth ms.
  (C0\=C1, (channel(C0, C1, Lat, Bwidth), depends_on(T0, T1, Data)) ->
      Com = Lat %+ Data/Bwidth
    ;
      Com = 0
  ),
  Tmp1 is Tmp + ET0 + Com,
  sumPCs([process_cost(T1, C1, ET1)|PCs], Tmp1, Res), !.

sumPCs([process_cost(T1, C1, ET1)|_], Tmp, Res):- Res is Tmp + ET1.
