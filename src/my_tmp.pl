%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_one
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



find_heuristically(Res):-
  find_dag(Res), !.

find_heuristically(Res):-
  find_batch_sol(Res), !.

find_dag(Res):-
  findall(Edge,
    make_dep_edge(Edge),
  EdgesArray),
  EdgesArray \= [],
  flatten(EdgesArray, Edges),
  vertices_edges_to_ugraph([], Edges, Graph),
  findall(C, core(C), Cores),
  length(Cores, LenCores),

  find_dag(Graph, LenCores, 0, [], Res).


find_dag([], _, _, ScheduleList, ScheduleList).

find_dag(Graph, LenCores, Counter, ScheduleList, Res):-
  get_core(Counter, LenCores, Core),
  NewCounter is Counter + 1,
  critical_path(Graph, MaxPath),
  add_path_to_schedule(MaxPath, Core, ScheduleList, NewScheduleList),
  del_vertices(Graph, MaxPath, NewGraph),
  find_dag(NewGraph, LenCores, NewCounter, NewScheduleList, Res).

get_core(Counter, LenCores, Core):-
  CoreNumberTmp is Counter mod LenCores, CoreNumber is CoreNumberTmp + 1,
  string_concat('c', CoreNumber, CoreStr), atom_string(Core, CoreStr).

find_batch_sol(Res):-
  findall(T, task(T), Tasks),
  et_sort(Tasks, SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], Res).

add_pc_list_to_schedule([], ScheduleList, ScheduleList).

add_pc_list_to_schedule([T|Ts], ScheduleList, Res):-
  findall(C, core(C), Cores),
  maplist(
    get_pc(T),
    Cores,
  PCs),
  minETList(ScheduleList, PCs, process_cost(Task, Core, MinTime)),
  add_pc_to_schedule(process_cost(Task, Core, MinTime), ScheduleList,
                     NewScheduleList),

  add_pc_list_to_schedule(Ts, NewScheduleList, Res).

get_pc( T, C, process_cost(T, C, ET)):- process_cost(T, C, ET).

% checks candidates for a solution. chooses the task and cpu which
% increases the total execution time the least
minETList(_, [X], X) :- !.
minETList(ScheduleList, [process_cost(T1, C1, ET1),
                         process_cost(T2, C2, ET2)|Tail], Res):-

  add_pc_to_schedule(process_cost(T1, C1, ET1), ScheduleList,
                       NewScheduleList1),
  execution_time(NewScheduleList1, TotalET1),
  add_pc_to_schedule(process_cost(T2, C2, ET2), ScheduleList,
                       NewScheduleList2),
  execution_time(NewScheduleList2, TotalET2),


  ( (TotalET1 < TotalET2) ->
      minETList(ScheduleList, [process_cost(T1, C1, ET1)|Tail], Res)

  ;
      minETList(ScheduleList, [process_cost(T2, C2, ET2)|Tail], Res)
  ), !.

% preprocessing: sorts task according to their execution time
et_sort(Tasks, Res):-
  maplist(
    get_min_pc,
  Tasks, PCs),
  flatten(PCs, FlattedPCs),
  predsort(compareTime, FlattedPCs, SortedPCs),
  maplist(extractTask, SortedPCs, Res).

extractTask(process_cost(T, _, _), T).
extractET(C, T, ET):- process_cost(T, C, ET).

% helper function for maplist above
get_min_pc(T, Res):-
  findall(process_cost(T, C, Ti), process_cost(T, C, Ti), PCs),
  minProcList(PCs, Res).

% gets the maximum processing cost in a list
%  TODO: rename
minProcList([X], X) :- !.
minProcList([process_cost(T1, C1, Ti1),
             process_cost(T2, C2, Ti2)|Tail],
            Res):-
    ( Ti1 > Ti2 ->
        minProcList([process_cost(T2, C2, Ti2)|Tail], Res)
    ;
        minProcList([process_cost(T1, C1, Ti1)|Tail], Res)
    ).

add_path_to_schedule([], _, ScheduleList, ScheduleList).

add_path_to_schedule([T|Ts], Core, ScheduleList, Res):-
  process_cost(T, Core, ET),
  add_pc_to_schedule(process_cost(T, Core, ET), ScheduleList, NewScheduleList),
  add_path_to_schedule(Ts, Core, NewScheduleList, Res).

% add_min(PCs, ScheduleList, Res):-
%   add_min(PCs, ScheduleList, [], Res).
%
% add_min([], _, Tmp, Tmp).
%
% add_min([PC|PCs], ScheduleList, Tmp, Res):-
%   add_pc_to_schedule(PC, ScheduleList, TmpScheduleList),
%   execution_time(TmpScheduleList, ET),
%   append([res(PC, ET)], Tmp, NewTmp),
%   add_min(PCs, ScheduleList, NewTmp, Res).

% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).


% valid schedule --> valid core and valid tasks
isSchedule(schedule(C, [T|Ts])):- core(C), isTaskSet(T, Ts).
isTaskSet(T, [To|Ts]):- task(T), isTaskSet(To, Ts).
isTaskSet(T, []):- task(T).


%  TODO: is Solution?
% isSolution(solution(Schedules)):-
%   findall(T, task(T), AllTasks),
%   length(AllTasks, Len),
%   isSolution(Schedules, 0, Num),
%   Len = Num.

% isSolution([schedule(C, Ts)|Schdls], Counter, Res):-
%   length(Ts, NumTs),
%   Counter1 is Counter + NumTs,
%   isSolution(Schdls, Counter1, Res).
%
% isSolution([], Counter, Counter).


find_all_sol(Sols, LenSols):-
  findall(Sol, find_sol(Sol), Sols),
  length(Sols, LenSols).


find_sol(Solution):-
  dep_sort_tasks(SortedTs),
  find_sol(SortedTs, [], Solution).

% create randomized solution
find_sol([T|Ts], ScheduleList, Solution):-
  process_cost(T, C, ET),
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).

% add process_cost(T, C, ET) to the corresponding schedule
add_pc_to_schedule(process_cost(Task, Core, Cost), Sol, NewSol):-
  maplist(add(process_cost(Task, Core, Cost)), Sol, NewSolTmp),
  (Sol = NewSolTmp ->
    append(NewSolTmp, [schedule(Core, [Task])], NewSol);
    append(NewSolTmp, [], NewSol)
  ).

add(process_cost(Task, Core, _), schedule(Core, Tasks),
  schedule(Core, NewTasks)):-
      append(Tasks, [Task], NewTasks),
      !.

add(process_cost(_, _, _), schedule(Core1, Tasks), schedule(Core1, Tasks)).

% sorting based on dependencies, i.e. topological sorting
dep_sort_PCs(Sorted):-
  findall(T, make_pc_dep_edge(T), Ts),
  (Ts = [] ->
    fail
  ;
    flatten(Ts, Fl),
    vertices_edges_to_ugraph([], Fl, L),
    top_sort(L, Sorted), !
  ).

dep_sort_tasks(Res):-
  dep_sort_PCs(SortedPCs),
  % TODO
  extract_list(SortedPCs, SortedTs),
  list_to_set(SortedTs, Res).

extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).

extract_list([], Tmp, Tmp).

extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).

% compare by time with respect to dependencies
compareTime(>, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1<ET2.
compareTime(<, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1>=ET2.

find_all_dep_edges(Res):-
  findall(Edge, make_dep_edge(Edge), Res).

find_all_dep_pcs(T, Res):-
  process_cost(T, C, ET),
  findall(PC, make_pc_dep_edge([process_cost(T, C, ET)-PC]), Res), !.

make_dep_edge([EnablingTask-DepTask]):- depends_on(DepTask, EnablingTask, _).

make_pc_dep_edge([process_cost(EnablingTask, C0, ET0)-process_cost(DepTask, C1, ET1)]):-
  !, process_cost(EnablingTask, C0, ET0),
  process_cost(DepTask, C1, ET1),
  depends_on(DepTask, EnablingTask, _).

make_pc_dep_edge(ScheduleList, [process_cost(Task, C1, ET1)-process_cost(DepTask, _, _)]):-
  depends_on(DepTask, Task, _),
  getPC(ScheduleList, Task, process_cost(Task, C1, ET1)).


critical_path(Graph, MaxPath):-
  critical_path(Graph, [], MaxPath).

critical_path([], Paths, MaxPath):-
  longest_path(Paths, MaxPath).

critical_path(Graph, OldPaths, MaxPath):-
  top_sort(Graph, [Root|_]),
  findall(Path,
    traverse(Graph, [[Root]], Path),
  Paths),

  flatten(Paths, Nodes),
  del_vertices(Graph, Nodes, NewGraph),
  append(Paths, OldPaths, NewPaths),
  critical_path(NewGraph, NewPaths, MaxPath).


get_edges(Path, Res):-
  get_edges(Path, [], Res).

get_edges([_], Res, Res).

get_edges([T1, T2|Ts], Tmp, Res):-
  append([T1-T2], Tmp, NewTmp),
  get_edges([T2|Ts], NewTmp, Res).


execution_time(ScheduleList, Res):-
  create_graph(ScheduleList, Graph),
  Graph \= [],
  max_pc_path(Graph, ScheduleList, [], tuple(_, Res)), !.

execution_time(ScheduleList, Res):-
  execution_time(ScheduleList, 0, Res).

execution_time([schedule(C, Ts)|Tail], Max, Res):-
  maplist(extractET(C), Ts, ETs),
  sumlist(ETs, TotalET),
  ( (TotalET > Max) ->
      NewMax is TotalET
  ;
      NewMax is Max
  ),
  execution_time(Tail, NewMax, Res).

execution_time([], Max, Max).

% depth first search
dfs_search(_, [Current|_], Target, Path):-
  % children(Graph, Current, Children),
  Current = Target,
  reverse(Current, Path).

dfs_search(Graph, [Current|Rest], Target, Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  dfs_search(Graph, NewAgenda, Target, Res).

% traverse graph to leaves
traverse(Graph, [Current|_], Path):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path), !.

traverse(Graph, [Current|Rest], Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  traverse(Graph, NewAgenda, Res), !.


children(Graph, [Node|RestOfPath], Children):-
  findall([Child,Node|RestOfPath],
    neighbour(Graph, Node, Child),
  Children).

% traverse graph to leaves
traverse_no_path(Graph, [Current|_], Path):-
  only_children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path).

traverse_no_path(Graph, [Current|Rest], Res):-
  only_children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  traverse_no_path(Graph, NewAgenda, Res).

only_children(Graph, Node, Children):-
  findall(Child,
    neighbour(Graph, Node, Child),
  Children).

% get process_costs for schedule
getPCWrapper(Schedules, Tasks, Res):-
  getPCWrapper(Schedules, Tasks, [], Res).

getPCWrapper(_, [], Tmp, Tmp).

getPCWrapper(Schedules, [T|Ts], Tmp, Res):-
  getPC(Schedules, T, PC),
  append([PC], Tmp, Tmp1),
  getPCWrapper(Schedules, Ts, Tmp1, Res).

getPC([Schedule|_], Task, Res):- getPC(Schedule, Task, Res).

getPC([_|Tail], Task, Res):- getPC(Tail, Task, Res).

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

sumPCs([process_cost(_, _, ET1)|_], Tmp, Res):- Res is Tmp + ET1.


make_io_pc_edges(ScheduleList, Schedule, Res):-
  make_io_pc_edges(ScheduleList, Schedule, [], Res).


make_io_pc_edges(_, schedule(_, []), Edges, EdgesSet):-
  list_to_set(Edges, EdgesSet).


make_io_pc_edges(Schedules, schedule(C, [T|Ts]), Edges0, Res):-
  process_cost(T, C, ET),
  findall(PC-process_cost(T, C, ET),
    make_pc_dep_edge(Schedules, [PC-process_cost(T, C, ET)]),
  Edges),
  append( Edges, Edges0, Edges1),
  make_io_pc_edges(Schedules, schedule(C, Ts), Edges1, Res), !.


% make_io_pc_edges(Schedules, schedule(C, [T1,T2|Ts]), Edges0, Res):-
%   process_cost(T1, C, ET1),
%   process_cost(T2, C, ET2),
%   append(
%     [process_cost(T1, C, ET1)-process_cost(T2, C, ET2)],
%   Edges0, Edges),
%   append( Edges0, Edges, Edges1),
%   make_io_pc_edges(Schedules, schedule(C, [T2|Ts]), Edges1, Res).

% make_io_pc_edges(Schedules, schedule(C, [T]), Edges0, Edges):-
%   process_cost(T, C, ET),
%   %TODO:check deps
%   findall(PC-process_cost(T, C, ET),
%     make_pc_dep_edge(Schedules, [PC-process_cost(T, C, ET)]),
%   DepEdges),
%   append(Edges0, DepEdges, Edges1),
%   make_io_pc_edges(Schedules, schedule(C, []), Edges1, Edges).
%
%
% make_io_pc_edges(Schedules, schedule(C, [T1,T2|Ts]), Edges0, Res):-
%   process_cost(T1, C, ET1),
%   process_cost(T2, C, ET2),
%   append(
%     [process_cost(T1, C, ET1)-process_cost(T2, C, ET2)],
%   Edges0, Edges1),
%   findall(PC-process_cost(T2, C, ET2),
%     make_pc_dep_edge(Schedules, [PC-process_cost(T2, C, ET2)]),
%   OtherEdges),
%   append( OtherEdges, Edges1, Edges2),
%   make_io_pc_edges(Schedules, schedule(C, [T2|Ts]), Edges2, Res).


make_io_pc_edges(_, schedule(_, [_]), Edges, Edges).

create_graph(ScheduleList, Res):-
  create_graph(ScheduleList, ScheduleList, [], Res).

create_graph(_, [], Edges, Graph):-
  vertices_edges_to_ugraph([], Edges, Graph).

create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  make_io_pc_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
  append(Edges0, InnerEdges, Edges2),
  create_graph(ScheduleList, Tail, Edges2, Res).


neighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).


longest_path([P|Ps], Max) :-
    longest_path(Ps, P, Max).

longest_path([], Max, Max).

longest_path([P1|Ps], P0, Max) :-
  length(P0, LenP0),
  length(P1, LenP1),

  (LenP1 > LenP0 ->
    MaxTmp = P1
  ;
    MaxTmp = P0
  ),
  longest_path(Ps, MaxTmp, Max).

find_all_paths(Graph, PC, Res):-
  find_all_paths(Graph, PC, [], Res).


find_all_paths([], _, Paths, Paths).

find_all_paths(Graph, PC, Paths, Res):-
  traverse(Graph, [ [PC] ], Path),
  append(Paths, [Path], NewPaths),
  del_vertices(Graph, Path, NewGraph),
  find_all_paths(NewGraph, PC, NewPaths, Res).

max_pc_path(_, [], Tmp1, Res):- maxPCs(Tmp1, Res).

max_pc_path(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
  process_cost(T, C, ET),
  findall(Path,
    traverse(Graph, [ [process_cost(T, C, ET)] ], Path),
  Paths),
  % find_all_paths(Graph, process_cost(T, C, ET), Paths),
  maxPCs(Paths, tuple(MaxPath, _)),
  append([MaxPath], Tmp, Tmp1),
  max_pc_path(Graph, Tail, Tmp1, Res).

