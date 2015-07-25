%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_one
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


find_optimal(Res):-
  findall(tuple(Sol, ET), find_sol(Sol, ET), SolList),
  min_ET(SolList, Res).

min_ET([Min],Min).

min_ET([tuple(S0, ET0), tuple(_, ET1)|T], Min):-
    ET0 =< ET1,
    min_ET([tuple(S0, ET0)|T], Min).

min_ET([tuple(_, ET0), tuple(S1, ET1)|T], Min):-
    ET0 > ET1,
    min_ET([tuple(S1, ET1)|T], Min).

find_sol(Solution, ET):-
  dep_sort_tasks(SortedTs),
  SortedTs \= [],
  find_sol(SortedTs, [], Solution),
  execution_time(Solution, ET).

% TODO: fix later
% find_sol(Solution, ET):-
%   et_sort(SortedTs),
%   SortedTs \= [],
%   find_sol(SortedTs, [], Solution),
%   execution_time(Solution, ET).

find_sol([T|Ts], ScheduleList, Solution):-
  process_cost(T, C, ET),
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).

find_heuristically(Res):-
  find_dag_sol(Res), !.

find_heuristically(Res):-
  find_batch_sol(Res), !.

find_dag_sol(Res):-
  findall(Edge,
    make_dep_edge(Edge),
  EdgesArray),
  EdgesArray \= [],
  flatten(EdgesArray, Edges),
  vertices_edges_to_ugraph([], Edges, Graph),
  findall(C, core(C), Cores),
  length(Cores, LenCores),

  find_dag_sol(Graph, LenCores, 0, [], Res).

find_dag_sol([], _, _, ScheduleList, ScheduleList).

find_dag_sol(Graph, LenCores, Counter, ScheduleList, Res):-
  get_core(Counter, LenCores, Core),
  NewCounter is Counter + 1,
  critical_path(Graph, MaxPath),
  add_path_to_schedule(MaxPath, Core, ScheduleList, NewScheduleList),
  del_vertices(Graph, MaxPath, NewGraph),
  find_dag_sol(NewGraph, LenCores, NewCounter, NewScheduleList, Res).

get_core(Counter, LenCores, Core):-
  CoreNumberTmp is Counter mod LenCores, CoreNumber is CoreNumberTmp + 1,
  string_concat('c', CoreNumber, CoreStr), atom_string(Core, CoreStr).

find_batch_sol(Res):-
  et_sort(SortedTasks),
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
et_sort(Res):-
  findall(T, task(T), Tasks),
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
  process_cost(EnablingTask, C0, ET0),
  process_cost(DepTask, C1, ET1),
  depends_on(DepTask, EnablingTask, _).

make_pc_dep_edge(ScheduleList, [process_cost(Task, C1, ET1)-process_cost(DepTask, _, _)]):-
  depends_on(DepTask, Task, _),
  getPC(ScheduleList, Task, process_cost(Task, C1, ET1)), !.

make_root_edge(ScheduleList, [process_cost(Task, _, _)-process_cost(DepTask, C, ET)]):-
  depends_on(DepTask, Task,  _),
  getPC(ScheduleList, DepTask, process_cost(DepTask, C, ET)).

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
  % TODO: fix
  getPC(ScheduleList, t1, Root),

  findall(Paths,
    bfs([ [Root] ], Graph, Paths),
    AllPaths
  ),
  maxPCs(AllPaths, tuple(_, Res)).

% execution_time(ScheduleList, Res):-
%   execution_time(ScheduleList, 0, Res).
%
% execution_time([schedule(C, Ts)|Tail], Max, Res):-
%   maplist(extractET(C), Ts, ETs),
%   sumlist(ETs, TotalET),
%   ( (TotalET > Max) ->
%       NewMax is TotalET
%   ;
%       NewMax is Max
%   ),
%   execution_time(Tail, NewMax, Res).
%
% execution_time([], Max, Max).

% shortest path
traverse(Graph, [Current|_], Path):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path), !.

traverse(Graph, [Current|Rest], Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  traverse(Graph, NewAgenda, Res).


children(Graph, [Node|RestOfPath], Children):-
  !, findall([Child,Node|RestOfPath],
    neighbour(Graph, Node, Child),
  Children), !.

% TODO: broken, get process_costs for schedule
getPC([Schedule|_], Task, Res):- getPC(Schedule, Task, Res).


getPC(schedule(C, Ts), Task, process_cost(Task, C, ET)):-
  member(Task, Ts),
  process_cost(Task, C, ET), !.

getPC([_|Schedules], Task, Res):- getPC(Schedules, Task, Res).

getPC([], _, _):- fail.


maxPCs([], tuple([], 0)).
maxPCs([[PCs]|Tail], Res):-
    sumPCs(PCs, Max),
    maxPCs(Tail, tuple(PCs, Max), Res).

maxPCs([], MaxTuple, MaxTuple).
maxPCs([[ PCs ]|Tail], tuple(PCs0, Max0), Res):-
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

make_io_edges(ScheduleList, Schedule, Res):-
  make_io_edges(ScheduleList, Schedule, [], Res).

make_io_edges(_, schedule(_, []), Edges, EdgesSet):-
  list_to_set(Edges, EdgesSet).

make_io_edges(Schedules, schedule(C, [T]), Edges0, Edges):-
  process_cost(T, C, ET),
  findall(process_cost(T, C, ET)-PC,
    make_root_edge(Schedules, [process_cost(T, C, ET)-PC]),
  DepEdges),
  append(Edges0, DepEdges, Edges1),
  make_io_edges(Schedules, schedule(C, []), Edges1, Edges), !.


make_io_edges(Schedules, schedule(C, [T1,T2|Ts]), Edges0, Res):-
  process_cost(T1, C, ET1),
  process_cost(T2, C, ET2),
  append(
    [process_cost(T1, C, ET1)-process_cost(T2, C, ET2)],
  Edges0, Edges1),
  findall(process_cost(T1, C, ET1)-PC,
    make_root_edge(Schedules, [process_cost(T1, C, ET1)-PC]),
  OtherEdges),
  append( OtherEdges, Edges1, Edges2),
  make_io_edges(Schedules, schedule(C, [T2|Ts]), Edges2, Res), !.


% make_io_pc_edges(ScheduleList, Schedule, Res):-
%   make_io_pc_edges(ScheduleList, Schedule, [], Res).
%
%
% make_io_pc_edges(_, schedule(_, []), Edges, EdgesSet):-
%   list_to_set(Edges, EdgesSet).
%
%
% make_io_pc_edges(Schedules, schedule(C, [T|Ts]), Edges0, Res):-
%   process_cost(T, C, ET),
%   findall(PC-process_cost(T, C, ET),
%     make_pc_dep_edge(Schedules, [PC-process_cost(T, C, ET)]),
%   Edges),
%   append( Edges, Edges0, Edges1),
%   make_io_pc_edges(Schedules, schedule(C, Ts), Edges1, Res), !.
%
%
% make_io_pc_edges(_, schedule(_, [_]), Edges, Edges).

create_graph(ScheduleList, Res):-
  create_graph(ScheduleList, ScheduleList, [], Res).

create_graph(_, [], Edges, Graph):-
  vertices_edges_to_ugraph([], Edges, Graph).

create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  % TODO: change back to pc
  make_io_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
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
  % findall(Path,
  %   traverse(Graph, [ [process_cost(T, C, ET)] ], Path),
  % Paths),
  findall(Paths,
    bfs([ [ process_cost(T, C, ET) ] ], Graph, Paths),
    AllPaths
  ),
  flatten(AllPaths, FlatPaths),
  % traverse_delete(Graph, process_cost(T, C, ET), Paths),
  maxPCs(FlatPaths, tuple(MaxPath, _)),
  append([MaxPath], Tmp, Tmp1),
  max_pc_path(Graph, Tail, Tmp1, Res).


empty_queue([]).
queue_head(S, Tail, [S|Tail]).
queue_last_list(Last, Q1, Q2):-
  reverse(Q1, RevQ1),
  reverse(Q2, [Last|RevQ1]).


% bfs(S, Graph, Path):-
%   empty_queue(Q1),
%   queue_head([S], Q1, Q2),
%   bfs1(Q2, Graph, Path).
%
% bfs1(Q, Graph, [Kids, S|Tail]):-
%   queue_head([S|Tail],_,Q),
%   children(Graph, S, Kids),
%   children \= [].
%
% bfs1(Q1, Graph, Solution):-
%   queue_head([S|Tail], Q2, Q1),
%   findall([Succ, S|Tail],
%   (neighbour(Graph, S, Succ), \+member(Succ,Tail)),
%   NewPaths),
%   queue_last_list(NewPaths,Q2,Q3),
%   bfs1(Q3, Graph, Solution).


bfs([Current|Rest], Graph, Path):-
  children(Graph, Current, Children),
  append(Rest,Children,NewAgenda),
  bfs(NewAgenda, Graph, Path).

bfs([Current|_], Graph, [Result]):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Result).

test(Rs):-
  execution_time([schedule(c1,[t1,t7]),schedule(c2,[t3,t5,t4,t6,t2])], Rs).
  % create_graph([schedule(c1, [t1, t2, t3, t4, t5, t6, t7])], G),
  % findall(R, bfs([ [ process_cost(t1, c1, 10) ] ], G, R), Rs).
