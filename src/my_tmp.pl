%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_one
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


find_optimal(Res):-
  dep_sort_tasks(DepSortedTs),
  ( (DepSortedTs \= []) ->
    findall(tuple(Sol, ET), find_sol_wrapper(DepSortedTs, Sol, ET), SolList)
  ;
    et_sort(EtSortedTs),
    findall(tuple(Sol, ET), find_sol_wrapper(EtSortedTs, Sol, ET), SolList)
  ),
  min_ET(SolList, Res), !.

min_ET([Min],Min).

min_ET([tuple(S0, ET0), tuple(_, ET1)|T], Min):-
    ET0 =< ET1,
    min_ET([tuple(S0, ET0)|T], Min).

min_ET([tuple(_, ET0), tuple(S1, ET1)|T], Min):-
    ET0 > ET1,
    min_ET([tuple(S1, ET1)|T], Min).

find_sol_wrapper(SortedTs, Solution, ET):-
  find_sol(SortedTs, [], Solution),
  execution_time(Solution, ET).

find_sol([T|Ts], ScheduleList, Solution):-
  ( (depends_on(_,_,_)) -> process_cost(T, C, ET); getPC(process_cost(T, C, ET))),
  % process_cost(T, C, ET),
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).

find_heuristically(Res):-
  find_dag_sol(Res).

find_heuristically(Res):-
  find_batch_sol(Res), !.

find_dag_sol(Res):-
  find_dag_sol(0, [], Res).

find_dag_sol(_, Paths, Paths):-
  flatten(Paths, Visited),
  list_to_set(Visited, Tasks),
  findall(T, task(T), AllTasks),
  length(Tasks, L0),
  length(AllTasks, L1),
  L0 = L1.

find_dag_sol(Counter, Paths, Res):-
  flatten(Paths, Visited),
  get_root(Visited, Root),
  my_dfs(Visited, Root, [Root], NewPath),
  append(Paths, [NewPath], NewPaths),
  NewCounter is Counter + 1,
  find_dag_sol(NewCounter, NewPaths, Res).

assign_cpu(Ps, Res):-
  assign_cpu(Ps, 0, [], Res).

assign_cpu([], _, ScheduleList, ScheduleList).

assign_cpu([P|Ps], Counter, ScheduleList, Res):-
  findall(C, core(C), Cores),
  length(Cores, LenCores),
  NewCounter is Counter + 1,
  get_core(Counter, LenCores, Core),
  add_path_to_schedule(P, Core, ScheduleList, NewScheduleList),
  assign_cpu(Ps, NewCounter, NewScheduleList, Res).

get_core(Counter, LenCores, Core):-
  CoreNumberTmp is Counter mod LenCores, CoreNumber is CoreNumberTmp + 1,
  string_concat('c', CoreNumber, CoreStr), atom_string(Core, CoreStr).


get_root(Root):-
  get_root([], Root).

get_root(Visited, Root):-
  findall(A-B, depends_on(B, A, _), Edges),
  vertices_edges_to_ugraph([], Edges, Graph),
  del_vertices(Graph, Visited, NewGraph),
  top_sort(NewGraph, [Root|_]).

get_root(Graph, _):-
  top_sort(Graph, [_|[Candid|_]]),
  del_vertices(Graph, [Candid], NewGraph),
  get_root(NewGraph, Candid).


get_inv_root(Root):-
  findall(A-B, depends_on(B, A, _), Edges),
  vertices_edges_to_ugraph([], Edges, Graph),
  transpose(Graph, TransGraph),
  top_sort(TransGraph, [Root|_]).


find_batch_sol(Res):-
  et_sort(SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], Res).

add_pc_list_to_schedule([], ScheduleList, ScheduleList).

add_pc_list_to_schedule([T|Ts], ScheduleList, Res):-
  findall(C, core(C), Cores),
  maplist(
    map_pc(T),
    Cores,
  PCs),
  minETList(ScheduleList, PCs, process_cost(Task, Core, MinTime)),
  add_pc_to_schedule(process_cost(Task, Core, MinTime), ScheduleList,
                     NewScheduleList),

  add_pc_list_to_schedule(Ts, NewScheduleList, Res).

map_pc( T, C, process_cost(T, C, ET)):- process_cost(T, C, ET).

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
  flatten(Ts, Fl),
  vertices_edges_to_ugraph([], Fl, L),
  top_sort(L, Sorted), !.

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


compareLen(>, A, B):- length(A, ALen), length(B, BLen), ALen<BLen.
compareLen(<, A, B):- length(A, ALen), length(B, BLen), ALen>=BLen.

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
  getPC(ScheduleList, DepTask, process_cost(DepTask, C, ET)), !.

critical_path(Root, MaxPath):-
  critical_path([], Root, [], MaxPath).

% critical_path(Visited, Root, Paths, Paths):-
%     my_dfs(Visited, Root, [], []).

critical_path(Visited, Root, OldPaths, MaxPath):-
  % dfs(Graph, [[Root]], Path),
  % my_dfs(Visited, Root, [], Path),
  my_dfs(Visited, Root, [Root], Path),
  append(Path, Visited, NewVisited),
  append([Path], OldPaths, NewPaths),
  critical_path(NewVisited, Root, NewPaths, MaxPath).


get_edges(Path, Res):-
  get_edges(Path, [], Res).

get_edges([_], Res, Res).

get_edges([T1, T2|Ts], Tmp, Res):-
  append([T1-T2], Tmp, NewTmp),
  get_edges([T2|Ts], NewTmp, Res).


execution_time(ScheduleList, Res):-
  !, depends_on(_, _, _),

  create_graph(ScheduleList, Graph),
  top_sort(Graph, [Root|_]),
  % getPC(ScheduleList, t1, Root),
  bfs(Graph, [[Root]], AllPaths),
  maxPCs(AllPaths, tuple(_, Res)).

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

% % shortest path
dfs(Graph, [Current|_], Path):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path), !.

dfs(Graph, [Current|Rest], Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  dfs(Graph, NewAgenda, Res), !.

bfs(Graph, [Current|_], Result):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Result).

bfs(Graph, [Current|Rest], Path):-
  children(Graph, Current, Children),
  append(Rest,Children,NewAgenda),
  bfs(Graph, NewAgenda, Path).


children(Graph, [Node|RestOfPath], Children):-
  findall([Child,Node|RestOfPath],
    neighbour(Graph, Node, Child),
  Children).

neighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).

my_dfs(Visited, Task, Path, R):-
  depends_on(DepTask, Task, _),
  not(member(DepTask, Visited)),
  append(Path, [DepTask], NewPath),
  my_dfs(Visited, DepTask, NewPath, R).

my_dfs(_, _, Path, Path):- !.



getPC(process_cost(T, c1, ET)):- process_cost(T, c1, ET).
getPC(process_cost(T, c2, ET)):- process_cost(T, c2, ET).
getPC(process_cost(T, c3, ET)):- process_cost(T, c3, ET).
getPC(process_cost(T, c4, ET)):- process_cost(T, c4, ET).
getPC(process_cost(T, C, ET)):- process_cost(T, C, ET), !.

% TODO: broken, get process_costs for schedule
getPC([Schedule|_], Task, Res):- getPC(Schedule, Task, Res).


getPC(schedule(C, Ts), Task, process_cost(Task, C, ET)):-
  member(Task, Ts),
  process_cost(Task, C, ET), !.

getPC([_|Schedules], Task, Res):- getPC(Schedules, Task, Res), !.

getPC([], _, _):- fail.


maxPCs([], tuple([], 0)):- !.
maxPCs([[PCs]|Tail], Res):-
    sumPCs(PCs, Max),
    maxPCs(Tail, tuple(PCs, Max), Res), !.

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
  findall(process_cost(T2, C, ET1)-PC,
    make_root_edge(Schedules, [process_cost(T2, C, ET1)-PC]),
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
  vertices_edges_to_ugraph([], Edges, Graph), !.

create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  % TODO: change back to pc
  make_io_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
  append(Edges0, InnerEdges, Edges2),
  create_graph(ScheduleList, Tail, Edges2, Res).



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
  dfs(Graph, [ [PC] ], Path),
  append(Paths, [Path], NewPaths),
  del_vertices(Graph, Path, NewGraph),
  find_all_paths(NewGraph, PC, NewPaths, Res).

% max_pc_path(_, [], Tmp1, Res):- maxPCs(Tmp1, Res).
%
% max_pc_path(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
%
%   process_cost(T, C, ET),
%   % findall(Path,
%   %   dfs(Graph, [ [process_cost(T, C, ET)] ], Path),
%   % Paths),
%   findall(Paths,
%     bfs([ [ process_cost(T, C, ET) ] ], Graph, Paths),
%     AllPaths
%   ),
%   flatten(AllPaths, FlatPaths),
%   % traverse_delete(Graph, process_cost(T, C, ET), Paths),
%   maxPCs(FlatPaths, tuple(MaxPath, _)),
%   append([MaxPath], Tmp, Tmp1),
%   max_pc_path(Graph, Tail, Tmp1, Res).


make_hg(process_cost(T, C, ET)-process_cost(TD, CD, ETD)):-
  process_cost(T, C, ET),
  process_cost(TD, CD, ETD),
  depends_on(TD, T, _).



map_task(process_cost(T, C, ET), T):- process_cost(T, C, ET).

delete_task_edges(Graph, [], Graph).
delete_task_edges(Graph, [T|Ts], Res):-
  findall(process_cost(T, C, ET), process_cost(T, C, ET), PCs),
  del_vertices(Graph, PCs, NewGraph),
  delete_task_edges(NewGraph, Ts, Res).

ex_init([], Max, Max).
ex_init(ScheduleList, R):-
  % getPC(ScheduleList, t7, RootPC),
  get_inv_root(Leave),

  findall(ET,
    ex(Leave, ScheduleList, ET),
    ETs
  ),
  maxList(ETs, R).

ex(T, ScheduleList, ET):-
  getPC(ScheduleList, T, process_cost(_, _, Cost)),
  findall(T0, depends_on(T, T0, _ ), Ts),
  Ts = [],
  ET is Cost.

ex(T, ScheduleList, ET):-
  getPC(ScheduleList, T, process_cost(_, _, Cost)),
  depends_on(T, T0, _ ),
  ex(T0, ScheduleList, ET0),
  ET is ET0 + Cost.

test(ET):-
  find_heuristically(List),
  predsort(compareLen, List, SList),
  assign_cpu(SList, ScheduleList),
  ex_init(ScheduleList, ET).
  %
  % top_sort(Graph, [Root|_]),
  % dfs(Graph, [[Root]], Path),
  % maplist(
  %   map_task,
  %   Path,
  % Tasks),
  % delete_task_edges(Graph, Tasks, NewGraph),
  % append(Paths, [Path], NewPaths),
  % test(NewGraph, NewPaths, Res).

  % find_heuristically(List),
  % assign_cpu(List, ScheduleList),
  % create_graph(ScheduleList, G),
  % bfs(G, [[process_cos(t1, c1, 10)]], R).
  % top_sort(G, [Root|_]),
  % findall(T1-T2, depends_on(T2, T1, _), Edges),
  % vertices_edges_to_ugraph([], Edges, Graph),
  % bfs(Graph, [[t1]], R).

  % execution_time(ScheduleList, ET).
  % create_graph([schedule(c1, [t1, t2, t3, t4, t5, t6, t7])], G),
  % findall(R, bfs([ [ process_cost(t1, c1, 10) ] ], G, R), Rs).
  %
%%
