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
  min_ET(SolList, tuple(Res, _)), !.

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
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).


find_heuristically(Res):-
  find_dag_sol(Res).

find_heuristically(Res):-
  find_batch_sol(Res), !.


find_dag_sol(Res):- find_dag_sol(0, [], Res).

find_dag_sol(_, List, ScheduleList):-
  flatten(List, Visited),
  list_to_set(Visited, Tasks),
  findall(T, task(T), AllTasks),
  length(Tasks, L0), length(AllTasks, L1),
  L0 = L1,
  predsort(compareLen, List, SList),
  assign_cpu(SList, ScheduleList).

find_dag_sol(Counter, Paths, Res):-
  flatten(Paths, Visited),
  get_root(Visited, Root),
  my_dfs(Visited, Root, [Root], NewPath),
  append(Paths, [NewPath], NewPaths),
  NewCounter is Counter + 1,
  find_dag_sol(NewCounter, NewPaths, Res).


assign_cpu(Ps, Res):-
  findall(schedule(C, []), core(C), ScheduleList),
  assign_cpu(Ps, 0, ScheduleList, Res).

assign_cpu([], _, ScheduleList, ScheduleList).

assign_cpu([P|Ps], Counter, ScheduleList, Res):-
  % get_valid_core_number(ScheduleList, P, Core),
  get_mod_core(Counter, Core),
  NewCounter is Counter + 1,
  create_graph(ScheduleList, G),
  top_sort(G, _),
  add_path_to_schedule(P, Core, ScheduleList, NewScheduleList),
  assign_cpu(Ps, NewCounter, NewScheduleList, Res).

assign_cpu([P|Ps], Counter, ScheduleList, Res):-
  NewCounter is Counter + 1,
  create_graph(ScheduleList, G),
  top_sort(G, _),
  add_path_to_schedule(P, c1, ScheduleList, NewScheduleList),
  assign_cpu(Ps, NewCounter, NewScheduleList, Res), !.


get_mod_core(Counter, Core):-
  findall(C, core(C), Cores), length(Cores, LenCores),
  CoreNumberTmp is Counter mod LenCores, CoreNumber is CoreNumberTmp + 1,
  string_concat('c', CoreNumber, CoreStr), atom_string(Core, CoreStr).

get_core_nr(Core, Nr):-
  atom_string(Core, CoreStr),
  string_concat('c', NrStr, CoreStr),
  atom_number(NrStr, Nr).

get_task_nr(Task, Nr):-
  atom_string(Task, TaskStr),
  string_concat('t', NrStr, TaskStr),
  atom_number(NrStr, Nr).


get_rndm_core(Cores, Core):-
  findall(C, core(C), Cores), length(Cores, LenCores),
  random(0, LenCores, CoreNumber),
  get_mod_core(CoreNumber, Core).


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
extractET(tuple(_, ET), ET).

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

add_path_to_schedule(P, Core, ScheduleList, R):-
  add_path_to_schedule(P, Core, ScheduleList, [], R).

add_path_to_schedule(P1, Core, [], Tmp, R):-
  append([schedule(Core, P1)], Tmp, TmpScheduleList),
  predsort(compare_core_nr, TmpScheduleList, R).

add_path_to_schedule(P1, Core, [schedule(Core, P0)], Tmp, R):-
  append(P0, P1, NewP),
  predsort(compare_task_nr, NewP, SortedP),
  append(Tmp, [schedule(Core, SortedP)], TmpScheduleList),
  predsort(compare_core_nr, TmpScheduleList, R).

add_path_to_schedule(P1, Core, [schedule(Core, P0)|Rest], Tmp, R):-
  append(P0, P1, NewP),
  predsort(compare_task_nr, NewP, SortedP),
  append(Tmp, [schedule(Core, SortedP)], Tmp1),
  append(Tmp1, Rest, TmpScheduleList),
  predsort(compare_core_nr, TmpScheduleList, R).

add_path_to_schedule(P1, C0, [schedule(C1, P0)|Rest], Tmp, R):-
  add_path_to_schedule(P1, C0, Rest, [schedule(C1, P0)|Tmp], R), !.


% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).


% valid schedule --> valid core and valid tasks
isSchedule(schedule(C, [T|Ts])):- core(C), isTaskSet(T, Ts).
isTaskSet(T, [To|Ts]):- task(T), isTaskSet(To, Ts).
isTaskSet(T, []):- task(T).

isSolution(solution(ScheduleList)):-
  create_graph(ScheduleList, Graph),
  top_sort(Graph,_).

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

% sorting based on no_deps, i.e. topological sorting
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

% compare by time with respect to no_deps
compareTime(>, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1<ET2.
compareTime(<, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1>=ET2.

compareLen(>, A, B):- length(A, ALen), length(B, BLen), ALen<BLen.
compareLen(<, A, B):- length(A, ALen), length(B, BLen), ALen>=BLen.

compareLast(<, A, B):-
  reverse(A, [LastA|_]),
  reverse(B, [LastB|_]),
  get_task_nr(LastA, LastANr),
  get_task_nr(LastB, LastBNr),
  LastANr<LastBNr.

compareLast(>, A, B):-
  reverse(A, [LastA|_]),
  reverse(B, [LastB|_]),
  get_task_nr(LastA, LastANr),
  get_task_nr(LastB, LastBNr),
  LastANr>=LastBNr.

compareLast1(<, [process_cost(FirstA, _, _)|_], [process_cost(FirstB, _, _)|_]):-
  get_task_nr(FirstA, FirstANr),
  get_task_nr(FirstB, FirstBNr),
  FirstANr<FirstBNr.

compareLast1(>, [process_cost(FirstA, _, _)|_], [process_cost(FirstB, _, _)|_]):-
  get_task_nr(FirstA, FirstANr),
  get_task_nr(FirstB, FirstBNr),
  FirstANr>=FirstBNr.


compare_task_nr(<, A, B):-
  get_task_nr(A, ANr),
  get_task_nr(B, BNr),
  ANr<BNr.

compare_task_nr(>, A, B):-
  get_task_nr(A, ANr),
  get_task_nr(B, BNr),
  ANr>=BNr.

compare_pc_task_nr(<, process_cost(A,_, _), process_cost(B, _, _)):-
  get_task_nr(A, ANr),
  get_task_nr(B, BNr),
  ANr<BNr.

compare_pc_task_nr(>, process_cost(A,_, _), process_cost(B, _, _)):-
  get_task_nr(A, ANr),
  get_task_nr(B, BNr),
  ANr>=BNr.

compare_core_nr(<, schedule(C0, _), schedule(C1, _)):-
  get_core_nr(C0, C0Nr),
  get_core_nr(C1, C1Nr),
  C0Nr<C1Nr.

compare_core_nr(>, schedule(C0, _), schedule(C1, _)):-
  get_core_nr(C0, C0Nr),
  get_core_nr(C1, C1Nr),
  C0Nr>=C1Nr.

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

critical_path(Root, MaxPath):-
  critical_path([], Root, [], MaxPath).

critical_path(Visited, _, Paths, Paths):-
  list_to_set(Visited, UniqVisted),
  length(UniqVisted, LenUniq),
  findall(T, task(T), Ts),
  length(Ts, LenTs),
  LenTs = LenUniq.

critical_path(Visited, Root, OldPaths, MaxPath):-
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
  depends_on(_, _, _),
  create_graph(ScheduleList, Graph),
  longest_path(Graph, Res), !.

execution_time(ScheduleList, Res):-
  execution_time(ScheduleList, 0, Res), !.

execution_time([schedule(C, Ts)|Tail], Max, Res):-
  maplist(extractET(C), Ts, ETs),
  sumlist(ETs, TotalET),
  ((TotalET > Max) -> NewMax is TotalET; NewMax is Max),
  execution_time(Tail, NewMax, Res), !.

execution_time([], Max, Max):- !.

% % shortest path
dfs([Current|_], Graph, Path):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Path).

dfs([Current|Rest], Graph, Res):-
  children(Graph, Current, Children),
  append(Children, Rest, NewAgenda),
  dfs(NewAgenda, Graph, Res).

bfs([Current|Rest], Graph, Path):-
  children(Graph, Current, Children),
  append(Rest,Children,NewAgenda), !,
  bfs(NewAgenda, Graph, Path).

bfs([Current|_], Graph, [Result]):-
  children(Graph, Current, Children),
  Children = [],
  reverse(Current, Result), !.


set_cpu(Children, Res):-
  findall(C, core(C), Cores),
  set_cpu(Children, Cores, [], Res).

set_cpu([], _, Tmp, Tmp).
set_cpu([T|Ts], [C|Cs], Tmp, Res):-
  % get_rndm_core(Cs, C),
  % delete(Cs, C, NewCs),
  process_cost(T, C, ET),
  set_cpu(Ts, Cs, [process_cost(T, C, ET)|Tmp], Res).

children(Graph, [Node|RestOfPath], SortedChildren):-
  findall([Child, Node|RestOfPath],
    neighbour(Graph, Node, Child),
  Children),
  predsort(compareLast1, Children, SortedChildren).

only_children(Graph, Node, SortedChildren):-
  findall(Child,
    neighbour(Graph, Node, Child),
  Children),
  predsort(compare_pc_task_nr, Children, SortedChildren).

neighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).

my_dfs(Visited, Task, Path, R):-
  depends_on(DepTask, Task, _),
  append(Path, [DepTask], NewPath),
  my_dfs(Visited, DepTask, NewPath, R), !.

% TODO: check!!!
my_dfs(Visited, _, Path, R):- subtract(Path, Visited, R).


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

create_graph(ScheduleList, Res):-
  create_graph(ScheduleList, ScheduleList, [], Res).


create_graph(_, [], Edges, Graph):-
  vertices_edges_to_ugraph([], Edges, Graph), !.

create_graph(ScheduleList, [schedule(_, [])|Tail], Edges0, Res):-
  create_graph(ScheduleList, Tail, Edges0, Res).

create_graph(ScheduleList, [schedule(C, [T|Ts])|Tail], Edges0, Res):-
  % TODO: change back to pc
  make_io_edges(ScheduleList, schedule(C, [T|Ts]), InnerEdges),
  append(Edges0, InnerEdges, Edges2),
  create_graph(ScheduleList, Tail, Edges2, Res).


get_dist(PC, [tuple(PC, W)], W):- !.
get_dist(PC, [tuple(PC, W)|_], W):- !.

get_dist(PC, [tuple(_, _)|Rest], R):-
  get_dist(PC, Rest, R).

get_dist(_, [], _):- !.

get_comm_cost(process_cost(CurT, CurC, _), process_cost(NbT, NbC, _), Comm):-
  (CurC\=NbC, (channel(CurC, NbC, Lat, Bwidth), depends_on(NbT, CurT, Data)) ->
      Comm = Lat + Data/Bwidth
    ;
      Comm = 0
  ).


add_dist(_, [], _, OldTuples, Tmp, R):-
  my_diff(Tmp, OldTuples, R), !.

add_dist(tuple(Cur, CurW), [process_cost(NbT, NbC, NbET)|Nbs], [tuple(process_cost(NbT, NbC, NbET), OldW)|_], OldTuples, Tmp, R):-
  get_comm_cost(Cur, process_cost(NbT, NbC, NbET), Comm),
  W1 is CurW + NbET + Comm,
  (W1 > OldW -> WNew is W1; WNew is OldW),
  append([tuple(process_cost(NbT, NbC, NbET), WNew)], Tmp, NewTmp),
  add_dist(tuple(Cur, CurW), Nbs, OldTuples, OldTuples, NewTmp, R), !.

add_dist(tuple(Cur, CurW), [Nb|Nbs], [tuple(_, _)|Tuples], OldTuples, Tmp, R):-
  add_dist(tuple(Cur, CurW), [Nb|Nbs], Tuples, OldTuples, Tmp, R).


init_length_to(process_cost(T, C, ET), tuple(process_cost(T, C, ET), ET)).
get_PC_tuple(tuple(PC, _), PC).


my_diff(List0, List1, R):-
  maplist(get_PC_tuple, List0, PCList0),
  maplist(get_PC_tuple, List1, PCList1),
  subtract(PCList1, PCList0, DiffList),
  copy_list(List1, DiffList, [], Rest),
  append(List0, Rest, R).


copy_list([], _, R, R):- !.

copy_list([tuple(PC, W)|Rest], DiffList, Tmp, R):-
  member(PC, DiffList),
  append([tuple(PC, W)], Tmp, NewTmp),
  copy_list(Rest, DiffList, NewTmp, R), !.

copy_list([_|Rest], DiffList, Tmp, R):-
  copy_list(Rest, DiffList, Tmp, R), !.


longest_path(G, Max):-
  top_sort(G, PCs),
  maplist(init_length_to, PCs, LengthTo),
  longest_path(G, PCs, LengthTo, List),
  maplist(extractET, List, ETList),
  max_list(ETList, Max).

longest_path(_, [], LengthTo, LengthTo):- !.

longest_path(G, [Cur|Rest], LenTo, R) :-
  get_dist(Cur, LenTo, CurW),
  only_children(G, Cur, Nbs),
  add_dist(tuple(Cur, CurW), Nbs, LenTo, LenTo, [], NewLenTo),
  longest_path(G, Rest, NewLenTo, R).

get_valid_core_number([schedule(C, [])], _, C).
get_valid_core_number([schedule(C, [])|_], _, C).

get_valid_core_number([schedule(C, Ts)|_], [First|_], C):-
  reverse(Ts, [Last|_]),
  get_task_nr(Last, LastNr),
  get_task_nr(First, FirstNr),
  FirstNr > LastNr.

get_valid_core_number([_|Rest], L, R):-
  get_valid_core_number(Rest, L, R).

hem(Res):-
  % findall(C, core(C), Cores),
  % init_scheduleList(Cores, [], [], TmpScheduleList),
  % findall(T0-T1, depends_on(T1, T0, _), Edges),
  % vertices_edges_to_ugraph([], Edges, G).
  dep_sort_tasks(SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], Res).

init_scheduleList([], ScheduleList, _, ScheduleList).
init_scheduleList([C|Cs], ScheduleList, Tmp, R):-
  task(T), not(member(T, Tmp)),
  add_path_to_schedule([T], C, ScheduleList, NewScheduleList),
  append([T], Tmp, NewTmp),
  init_scheduleList(Cs, NewScheduleList, NewTmp, R).


test(ET):-
  % hem(R).
  find_heuristically(ScheduleList),
  execution_time(ScheduleList, ET).
