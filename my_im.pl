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

delete_node(PC, [Path|Tail], Res):-
  delete_node(PC, [Path|Tail], [], Res).

delete_node(PC, [], Paths, Paths).

delete_node(PC, [Path|Tail], Paths0, Res):-
  del(PC, Path, NewPath),
  append(Paths0, [NewPath], Paths1),
  delete_node(PC, Tail, Paths1, Res).

% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).

minList([result(L, ET)|Ls], Min) :-
    minList(Ls, result(L, ET), Min).

minList([], Min, Min).
minList([result(L, ET)|Ls], result(L0, Min0), Min) :-
  (ET < Min0 ->
    Min1 = result(L, ET)
  ;
    Min1 = result(L0, Min0)
  ),
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

add1(PC, List, NewList):-
  append(PC, [List], NewList).

% get list with random element (not in TakenElems)
get_rndm_elem(List, RndmElem):-
  length(List, Len),
  Len1 is Len - 1,
  random_between(0, Len1, RndmIndex),
  nth0(RndmIndex, List, RndmElem), !.

%get_pc(ScheduleList, process_cost(Task, Core, MinTime)):-
  %% TODO: trick the cuts
  %findall(process_cost(Task, C, Time),
          %process_cost(Task, C, Time), ProcCosts),
  %minETList(ScheduleList, ProcCosts, process_cost(Task, Core, MinTime)),
  %% important, only one solution here.
  %process_cost(Task, Core, MinTime), !.
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

% required predicate to find the exact solution
find_optimal(result(Sol, ET)):-
  depSort(TaskList),
  find_one(TaskList, [], result(Sol, ET)), !.

find_one([], ScheduleList, result(ScheduleList, ET)):-
  execution_time(ScheduleList, ET), !.

find_one([Task|Tasks], ScheduleList, result(ResSchedule, ET)):-
  % binds only, if stated more solutions possible (?)
  findall(process_cost(Task, C, Time),
    get_pc(process_cost(Task, C, Time)),
  ProcCosts),
  minETList(ScheduleList, ProcCosts, process_cost(Task, Core, MinTime)),
  add_to_schedule_list(process_cost(Task, Core, MinTime),
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
  findall(process_cost(T, C, Ti), get_pc(process_cost(T, C, Ti)), Tasks),
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

extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).
extract_list([], Tmp, Tmp).
extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).


make_edge(ScheduleList,
  [process_cost(Task, C1, ET1)-process_cost(DepTask, C2, ET2)]):-
  depends_on(DepTask, Task, _),
  getPC(ScheduleList, Task, process_cost(Task, C1, ET1)).

make_edge([EnablingTask-DepTask]):- depends_on(DepTask, EnablingTask, _).

% compare by time with respect to dependencies
compareTime(>, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1<ET2.

compareTime(<, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1>=ET2.

execution_time(ScheduleList, Res):-
  create_graph(ScheduleList, Graph),
  %findall(Path,
    %search_df(Graph, [[process_cost(start, null, null)]], Path),
    %Paths),
  %delete_node(process_cost(start, null, null), Paths, Paths1),
  %maxPCs(Paths1, tuple(_, Res)), !.
  search_df_wrapper(Graph, ScheduleList, [], Res), !.

% depth first search
search_df(Graph, [Current|_], Path):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path).

search_df(Graph, [Current|Rest], Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  search_df(Graph, NewAgenda, Res).

children(Graph, [Node|RestOfPath], Children):-
  findall([Child,Node|RestOfPath],
    isNeighbour(Graph, Node,Child),
  Children).


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


maxPCs([PCs|Tail], Res):-
    sumPCs(PCs, Max),
    maxPCs(Tail, tuple(PCs, Max), Res).

maxPCs([], MaxTuple, MaxTuple).
maxPCs([PCs|Tail], tuple(PCs0, Max0), Res):-
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
  (C0\=C1, (channel(C0, C1, Lat, Bwidth), depends_on(T1, T0, Data)) ->
      Com = Lat + Data/Bwidth
    ;
      Com = 0
  ),
  Tmp1 is Tmp + ET0 + Com,
  sumPCs([process_cost(T1, C1, ET1)|PCs], Tmp1, Res), !.

sumPCs([process_cost(T1, C1, ET1)|_], Tmp, Res):- Res is Tmp + ET1.

make_inner_edges(ScheduleList, Schedule, Res):-
  make_inner_edges(ScheduleList, Schedule, [], Res).

make_inner_edges(_, schedule(C, []), Edges, EdgesSet):-
  list_to_set(Edges, EdgesSet).
make_inner_edges(Schedules, schedule(C, [T]), Edges0, Edges):-
  process_cost(T, C, ET),
  %TODO:check deps
  findall(PC-process_cost(T, C, ET),
    make_edge(Schedules, [PC-process_cost(T, C, ET)]),
  DepEdges),
  append(Edges0, DepEdges, Edges1),
  make_inner_edges(Schedules, schedule(C, []), Edges1, Edges).


make_inner_edges(Schedules, schedule(C, [T1,T2|Ts]), Edges0, Res):-
  process_cost(T1, C, ET1),
  process_cost(T2, C, ET2),
  append(
    [process_cost(T1, C, ET1)-process_cost(T2, C, ET2)],
  Edges0, Edges1),
  findall(PC-process_cost(T2, C, ET2),
    make_edge(Schedules, [PC-process_cost(T2, C, ET2)]),
  OtherEdges),
  append( OtherEdges, Edges1, Edges2),
  make_inner_edges(Schedules, schedule(C, [T2|Ts]), Edges2, Res).


create_graph(ScheduleList, Res):-
  create_graph(ScheduleList, ScheduleList, [], Res).

create_graph(_, [], Edges, Graph):-
  vertices_edges_to_ugraph([], Edges, Graph).
create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  %process_cost(T, C, ET),
  %append(
    %[process_cost(start, null, null)-process_cost(T, C, ET)],
  %Edges0, Edges1),

  %append(Edges1, InnerEdges, Edges2),

  make_inner_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
  append(Edges0, InnerEdges, Edges2),
  create_graph(ScheduleList, Tail, Edges2, Res).

isNeighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).


find_my_optimal(MinET):-
  depSort(TaskList),
  findall(Sol,
    find_my_optimal(TaskList, [], Sol),
  AllSols),
  minList(AllSols, MinET).

find_my_optimal([], ScheduleList, result(ScheduleList, ET)):-
  execution_time(ScheduleList, ET).

find_my_optimal([Task|Tasks], ScheduleList, Res):-
  % binds only, if stated more solutions possible (?)
  find_optimal_pc(Task, ScheduleList, NewScheduleList),
  find_my_optimal(Tasks, NewScheduleList, Res).
  %isSolution(isSolution(C)).
  % TODO: get the schedule with the minimum execution time

find_optimal_pc(Task, ScheduleList, NewSol):-
  % crucial point, more solutions are possible for the following line
  get_pc(process_cost(Task, Core, Time)),
  process_cost(Task, Core, Time),
  add_to_schedule_list(process_cost(Task, Core, Time), ScheduleList,
  NewSol).

search_df_wrapper(Graph, [], Tmp1, Res):- maxPCs(Tmp1, tuple(_, Res)).

search_df_wrapper(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
    isNeighbour(Graph, Predecessor, process_cost(T, C, ET)),
    search_df_wrapper(Graph, Tail, Tmp, Res).

search_df_wrapper(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
  process_cost(T, C, ET),
  findall(Path,
    search_df(Graph, [[process_cost(T, C, ET)]], Path),
    Paths),
    maxPCs(Paths, tuple(MaxPath, _)),
    append([MaxPath], Tmp, Tmp1),
    search_df_wrapper(Graph, Tail, Tmp1, Res).


