%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% quite sketchy, but first important functions:
% find_heuristically, find_one
% check below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

sarkar(Res):-
  dep_sort_tasks([T|Ts]),
  process_cost(T, C, ET),
  add_to_schedule_list(process_cost(T, C, ET), [], TmpScheduleList),
  sarkar(Ts, TmpScheduleList, Res).


sarkar([], ScheduleList, Res).

sarkar([T|Ts], ScheduleList, Res):-
  find_all_dep_pcs(T, PCs),
  predsort(compareTime, PCs, SortedPCs),
  add_min(SortedPCs, ScheduleList, Res).
  % sarkar(Ts, ScheduleList, Res), !.


add_min(PCs, ScheduleList, Res):-
  add_min(PCs, ScheduleList, [], Res).

add_min([], _, Tmp, Tmp).

add_min([PC|PCs], ScheduleList, Tmp, Res):-
  add_to_schedule_list(PC, ScheduleList, TmpScheduleList),
  execution_time(TmpScheduleList, ET),
  append([res(PC, ET)], Tmp, NewTmp),
  add_min(PCs, ScheduleList, NewTmp, Res).

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


find_all_sol(Sols, LenSols):-
  findall(Sol, find_sol(Sol), Sols),
  length(Sols, LenSols).


find_sol(Solution):-
  dep_sort_tasks(SortedTs),
  find_sol(SortedTs, [], Solution).

% create randomized solution
find_sol([T|Ts], ScheduleList, Solution):-
  process_cost(T, C, ET),
  add_to_schedule_list(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).

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

add(process_cost(_, _, _), schedule(Core1, Tasks), schedule(Core1, Tasks)).

% sorting based on dependencies, i.e. topological sorting
dep_sort_PCs(Sorted):-
  findall(T, make_dep_edge(T), Ts),
  (Ts = [] ->
    fail
  ;
    flatten(Ts, Fl),
    vertices_edges_to_ugraph([], Fl, L),
    top_sort(L, Sorted), !
  ).

dep_sort_tasks(Res):-
  dep_sort_PCs(SortedPCs),
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
compareTime(>, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1<ET2.
compareTime(<, process_cost(T1, _, ET1), process_cost(T2, _, ET2)):-
  ET1>=ET2.

find_all_dep_edges(Res):-
  findall(Edge, make_dep_edge(Edge), Res).


find_all_dep_pcs(T, Res):-
  process_cost(T, C, ET),
  findall(PC, make_dep_edge([process_cost(T, C, ET)-PC]), Res), !.


make_dep_edge([process_cost(EnablingTask, C0, ET0)-process_cost(DepTask, C1, ET1)]):-
  !, process_cost(EnablingTask, C0, ET0),
  process_cost(DepTask, C1, ET1),
  depends_on(DepTask, EnablingTask, _).

make_dep_edge(ScheduleList, [process_cost(Task, C1, ET1)-process_cost(DepTask, C2, ET2)]):-
  depends_on(DepTask, Task, _),
  getPC(ScheduleList, Task, process_cost(Task, C1, ET1)).

make_dep_edge([EnablingTask-DepTask]):- depends_on(DepTask, EnablingTask, _).

make_dep_edge_rev([process_cost(DepTask, C1, ET1)-process_cost(EnablingTask, C0, ET0)]):-
  !, process_cost(EnablingTask, C0, ET0),
  process_cost(DepTask, C1, ET1),
  depends_on(DepTask, EnablingTask, _).


critical_path(MaxPath):-
  process_cost(t1, C1, ET1),
  findall(Edge,
    make_dep_edge(Edge),
  EdgesArray),
  flatten(EdgesArray, Edges),
  vertices_edges_to_ugraph([], Edges, Graph),
  % children(Graph, [process_cost(t1, C1, ET1)], Res).
  findall(Path,
    search_df(Graph, [[process_cost(t1, C1, ET1)]], Path),
    Paths),
  maxPCs(Paths, tuple(MaxPath, _)).

execution_time(ScheduleList, Res):-
  create_graph(ScheduleList, Graph),
  max_path(Graph, ScheduleList, [], tuple(_, Res)), !.

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

getPC([Schedule|Tail], Task, Res):- getPC(Schedule, Task, Res).

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

sumPCs([process_cost(T1, C1, ET1)|_], Tmp, Res):- Res is Tmp + ET1.

make_io_edges(ScheduleList, Schedule, Res):-
  make_io_edges(ScheduleList, Schedule, [], Res).

make_io_edges(_, schedule(C, []), Edges, EdgesSet):-
  list_to_set(Edges, EdgesSet).

make_io_edges(Schedules, schedule(C, [T]), Edges0, Edges):-
  process_cost(T, C, ET),
  %TODO:check deps
  findall(PC-process_cost(T, C, ET),
    make_dep_edge(Schedules, [PC-process_cost(T, C, ET)]),
  DepEdges),
  append(Edges0, DepEdges, Edges1),
  make_io_edges(Schedules, schedule(C, []), Edges1, Edges).


make_io_edges(Schedules, schedule(C, [T1,T2|Ts]), Edges0, Res):-
  process_cost(T1, C, ET1),
  process_cost(T2, C, ET2),
  append(
    [process_cost(T1, C, ET1)-process_cost(T2, C, ET2)],
  Edges0, Edges1),
  findall(PC-process_cost(T2, C, ET2),
    make_dep_edge(Schedules, [PC-process_cost(T2, C, ET2)]),
  OtherEdges),
  append( OtherEdges, Edges1, Edges2),
  make_io_edges(Schedules, schedule(C, [T2|Ts]), Edges2, Res).


create_graph(ScheduleList, Res):-
  create_graph(ScheduleList, ScheduleList, [], Res).

create_graph(_, [], Edges, Graph):-
  vertices_edges_to_ugraph([], Edges, Graph).

create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  make_io_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
  append(Edges0, InnerEdges, Edges2),
  create_graph(ScheduleList, Tail, Edges2, Res).


isNeighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).


max_path(Graph, [], Tmp1, Res):- maxPCs(Tmp1, Res).

max_path(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
    isNeighbour(Graph, Predecessor, process_cost(T, C, ET)),
    max_path(Graph, Tail, Tmp, Res).

max_path(Graph, [schedule(C, [T|_])|Tail], Tmp, Res):-
  process_cost(T, C, ET),
  findall(Path,
    search_df(Graph, [[process_cost(T, C, ET)]], Path),
    Paths),
    maxPCs(Paths, tuple(MaxPath, _)),
    append([MaxPath], Tmp, Tmp1),
    max_path(Graph, Tail, Tmp1, Res).

