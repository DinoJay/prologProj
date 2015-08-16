%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%  My Implementation %%%
% find_optimal
% generates all possible solutions and picks the best
%
% find_heuristically
% based on "chooseMinPC", selecting the minimal processing_cost() which increases
% the ET the least and then adding it to the schedule with "add_pc_to_schedule"
%
% execution_time
% computes the ET with the longest path algorithm. It uses caching to keep the
% ET of predecessors of the current task. execution_time traverses the schedule
% in topological ordering if possible. In each iteration the direct predecessors
% are looked up and updated with the ET of the current task.
% find_heuristically reuses this cache in its iterations,
% It constructs a solution by adding tasks one by one in topological sorting to
% schedule.
%
% Important helpers:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

pretty_print([]).
pretty_print([schedule(C, Ts)|R]):-
  write(C), write(": "),
  writeln(Ts),
  pretty_print(R).

speedup(solution(S), Speedup):-
  speedup(S, Speedup).

speedup(S, Speedup):-
  find_optimal_seq(SeqOp),
  execution_time(SeqOp, ETSeqOp),
  execution_time(S, SET),
  Speedup is ETSeqOp / SET.


% find optimal ET if one put all tasks on single cores
find_optimal_seq(Res):-
  dep_sort_tasks(DepSortedTs),
  findall(C, core(C), Cs),
  maplist(assign_core(DepSortedTs), Cs, SolList),
  min_ET(SolList, tuple(Schedule, _)),
  predsort(compare_core_nr, Schedule, Res), !.

find_optimal_seq(Res):-
  et_sort(SortedTs),
  findall(C, core(C), Cs),
  maplist(assign_core(SortedTs), Cs, SolList),
  min_ET(SolList, tuple(Res, _)), !.


assign_core(Ts, C, tuple([schedule(C, Ts)], ET)):- execution_time([schedule(C, Ts)], ET).


% find optimal solution for sets with dependencies
find_optimal(Res):-
  dep_sort_tasks(DepSortedTs),
  DepSortedTs \= [],
  findall(tuple(Sol, ET), find_sol_exec_tuple(DepSortedTs, Sol, ET), SolList),
  min_ET(SolList, tuple(Res, _)), !.


% find optimal solution for the two first batch sets
find_optimal(Res):-
  et_sort(EtSortedTs),
  findall(tuple(Sol, ET), find_sol_exec_tuple(EtSortedTs, Sol, ET), SolList),
  min_ET(SolList, tuple(Res, _)), !.


% find min schedule with minimum ET in a tuple list consisting of schedule and
% ET
min_ET([Min],Min).

min_ET([tuple(S0, ET0), tuple(_, ET1)|T], Min):-
    ET0 =< ET1,
    min_ET([tuple(S0, ET0)|T], Min).

min_ET([tuple(_, ET0), tuple(S1, ET1)|T], Min):-
    ET0 > ET1,
    min_ET([tuple(S1, ET1)|T], Min).

% helper predicate to constructuct one solution, used to above to find all
% possible solutions
find_sol_exec_tuple(SortedTs, Solution, ET):-
  find_sol(SortedTs, [], Solution),
  execution_time(Solution, ET).

find_sol([T|Ts], ScheduleList, Solution):-
  ( (depends_on(_,_,_)) -> process_cost(T, C, ET); getPC(process_cost(T, C, ET))),
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, TmpScheduleList),
  find_sol(Ts, TmpScheduleList, Solution).

find_sol([], ScheduleList, ScheduleList).


% find heuristically for batch sets, sorts task by ET
find_heuristically(Res):-
  et_sort(SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], [], Schedule),
  predsort(compare_core_nr, Schedule, Res), !.

% dependency version
find_heuristically(Res):-
  depends_on(_, _, _),
  dep_sort_tasks(SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], [], Schedule),
  predsort(compare_core_nr, Schedule, Res), !.

% add task one by one to the schedule governed by least impact metric.
% Meaning, select the task which increases the current schedule the least
% efficiency is preserved by keeping the LenTo list which holds the ET for all
% predecessors
add_pc_list_to_schedule([], _, ScheduleList, ScheduleList).

add_pc_list_to_schedule([T|Ts], LenTo, ScheduleList, Res):-
  findall(C, core(C), Cores),
  maplist(map_pc(T), Cores, PCs),

  chooseMinPC(ScheduleList, LenTo, PCs, NewLenTo, process_cost(Task, Core, MinTime)),
  add_pc_to_schedule(process_cost(Task, Core, MinTime), ScheduleList, NewScheduleList),
  % important, if it fails jump to the batch version below
  top_sorted(NewScheduleList),
  add_pc_list_to_schedule(Ts, NewLenTo, NewScheduleList, Res).

map_pc( T, C, process_cost(T, C, ET)):- process_cost(T, C, ET).

% batch version
% checks candidates for a solution. chooses the task and cpu which
% increases the total execution time the least
chooseMinPC(ScheduleList, LenTo, [process_cost(T, C, ET)], NewLenTo, process_cost(T, C, ET)):-
  append(LenTo, [tuple(process_cost(T, C, ET), ET)], LenTo1),
  add_pc_to_schedule(process_cost(T, C, ET), ScheduleList, NewScheduleList1),
  incr_exec_time(NewScheduleList1, process_cost(T, C, ET), LenTo1, NewLenTo, _).

chooseMinPC(ScheduleList, LenTo, [process_cost(T0, C0, ET0), process_cost(T1, C1, ET1)|Tail], ResLenT0, Res):-

  append(LenTo, [tuple(process_cost(T0, C0, ET0), ET0)], LenTo0),
  add_pc_to_schedule(process_cost(T0, C0, ET0), ScheduleList, NewScheduleList0),
  incr_exec_time(NewScheduleList0, process_cost(T0, C0, ET0), LenTo0, _, TotalET0),

  append(LenTo, [tuple(process_cost(T1, C1, ET1), ET1)], LenTo1),
  add_pc_to_schedule(process_cost(T1, C1, ET1), ScheduleList, NewScheduleList1),
  incr_exec_time(NewScheduleList1, process_cost(T1, C1, ET1), LenTo1, _, TotalET1),
  ( (TotalET0 < TotalET1) ->
      chooseMinPC(ScheduleList, LenTo, [process_cost(T0, C0, ET0)|Tail], ResLenT0, Res)

  ;
      chooseMinPC(ScheduleList, LenTo, [process_cost(T1, C1, ET1)|Tail], ResLenT0, Res)
  ), !.

% sorts task according to their execution time
et_sort(Res):-
  findall(T, task(T), Tasks),
  maplist(
    get_min_pc,
  Tasks, PCs),
  flatten(PCs, FlattedPCs),
  predsort(compareTime, FlattedPCs, SortedPCs),
  maplist(extractTask, SortedPCs, Res).

% helper
extractTask(process_cost(T, _, _), T).
extractET(C, T, ET):- process_cost(T, C, ET).
extractET(tuple(_, ET), ET).

% helper function for maplist above
get_min_pc(T, Res):-
  findall(process_cost(T, C, Ti), process_cost(T, C, Ti), PCs),
  minPCList(PCs, Res).

% gets the minimum processing cost in a list
minPCList([X], X) :- !.
minPCList([process_cost(T1, C1, Ti1),
             process_cost(T2, C2, Ti2)|Tail],
            Res):-
    ( Ti1 > Ti2 ->
        minPCList([process_cost(T2, C2, Ti2)|Tail], Res)
    ;
        minPCList([process_cost(T1, C1, Ti1)|Tail], Res)
    ).

% find max in list
maxList([A],A).
maxList([A|List],Max):- maxList(List,Max1),
                        (A>=Max1, Max=A; A<Max1, Max=Max1).

% checks if a solution is valid, applies topological sort, occurency check that
% every core is only listed once and that all tasks are assigned
isSolution(solution(ScheduleList)):-
  depends_on(_, _, _),
  !,
  create_graph(ScheduleList, Graph),
  top_sort(Graph, Vertices),
  length(Vertices, LenVertices),
  findall(T, task(T), Ts),
  length(Ts, LenTs),
  LenVertices = LenTs,
  maplist(extract_core, ScheduleList, Cores),
  findall(C, core(C), AllCores),
  contains(AllCores, Cores), !.

% differes slightly for the batch sets, no top_sort!
isSolution(solution(ScheduleList)):-
  maplist(extract_tasks, ScheduleList, UnflatTs),
  flatten(UnflatTs, Ts),
  length(Ts, LenTs),
  findall(T, task(T), AllTs),
  length(AllTs, LenAllTs),
  LenTs = LenAllTs,
  maplist(extract_core, ScheduleList, Cores),
  findall(C, core(C), AllCores),
  contains(AllCores, Cores), !.

% List contains element?
contains(_, []):- !.
contains(List, [H|T]):- member(H, List), contains(List, T), !.

extract_core(schedule(C, _), C).
extract_tasks(schedule(_, Ts), Ts).

% checks if schedule is topologicaly sorted
top_sorted(ScheduleList):-
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
  extract_list(SortedPCs, SortedTs),
  list_to_set(SortedTs, Res).

% extract tasks from processing_cost list
extract_list(List, Res):-
  extract_list(List, [], OrderedList),
  reverse(OrderedList, Res).

extract_list([], Tmp, Tmp).

extract_list([process_cost(T, _, _)|PCs], Tmp, Res):-
  append([T], Tmp, Tmp1),
  extract_list(PCs, Tmp1, Res).

% get integer number from core
get_core_nr(Core, Nr):-
  atom_string(Core, CoreStr),
  string_concat('c', NrStr, CoreStr),
  atom_number(NrStr, Nr).

compare_core_nr(<, schedule(C0, _), schedule(C1, _)):-
  get_core_nr(C0, C0Nr),
  get_core_nr(C1, C1Nr),
  C0Nr<C1Nr.

compare_core_nr(>, schedule(C0, _), schedule(C1, _)):-
  get_core_nr(C0, C0Nr),
  get_core_nr(C1, C1Nr),
  C0Nr>=C1Nr.

% compare by time with respect to no_deps
compareTime(>, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1<ET2.
compareTime(<, process_cost(_, _, ET1), process_cost(_, _, ET2)):-
  ET1>=ET2.

compareLen(>, A, B):- length(A, ALen), length(B, BLen), ALen<BLen.
compareLen(<, A, B):- length(A, ALen), length(B, BLen), ALen>=BLen.

% helper functions to create edges from dependencies
make_pc_dep_edge([process_cost(EnablingTask, C0, ET0)-process_cost(DepTask, C1, ET1)]):-
  depends_on(DepTask, EnablingTask, _),
  process_cost(EnablingTask, C0, ET0),
  process_cost(DepTask, C1, ET1).

make_pc_dep_edge(ScheduleList, [process_cost(Task, C1, ET1)-process_cost(DepTask, _, _)]):-
  depends_on(DepTask, Task, _),
  getPC(ScheduleList, Task, process_cost(Task, C1, ET1)), !.

make_edge(ScheduleList, [process_cost(Task, _, _)-process_cost(DepTask, C, ET)]):-
  depends_on(DepTask, Task,  _),
  getPC(ScheduleList, DepTask, process_cost(DepTask, C, ET)).


% execution time with the use of temporary solutions
incr_exec_time(ScheduleList, PC, LenTo, NewLenTo, Res):-
  depends_on(_, _, _),
  create_graph(ScheduleList, Graph),
  predecessors(PC, Graph, PCs),
  my_longest_path(Graph, PCs, LenTo, NewLenTo, Res), !.

incr_exec_time(ScheduleList, _, _, _, Res):-
  simple_exec_time(ScheduleList, Res).

predecessors(Node, Graph, Predecessors):-
  transpose(Graph, RevGraph),
  neighbours(Node, RevGraph, Predecessors).

predecessors(_, _, []).

% execution time
execution_time(solution(ScheduleList), Res):-
  execution_time(ScheduleList, Res).

execution_time(ScheduleList, Res):-
  create_graph(ScheduleList, Graph),
  top_sort(Graph, _),
  longest_path(Graph, Res), !.

% faster execution time for batch sets
simple_exec_time(ScheduleList, Res):-
  simple_exec_time(ScheduleList, 0, Res), !.

simple_exec_time([schedule(C, Ts)|Tail], Max, Res):-
  maplist(extractET(C), Ts, ETs),
  sumlist(ETs, TotalET),
  ((TotalET > Max) -> NewMax is TotalET; NewMax is Max),
  simple_exec_time(Tail, NewMax, Res), !.

simple_exec_time([], Max, Max):- !.

children(Graph, Node, Children):-
  findall(Child,
    neighbour(Graph, Node, Child),
  Children).

neighbour(Graph, Node, Neighbour):-
  neighbours(Node, Graph, AllNeighbours),
  member(Neighbour, AllNeighbours).

% circumvent cut which is in some sets and restricts backtracking
getPC(process_cost(T, c1, ET)):- process_cost(T, c1, ET).
getPC(process_cost(T, c2, ET)):- process_cost(T, c2, ET).
getPC(process_cost(T, c3, ET)):- process_cost(T, c3, ET).
getPC(process_cost(T, c4, ET)):- process_cost(T, c4, ET).
getPC(process_cost(T, C, ET)):- process_cost(T, C, ET), !.

getPC([Schedule|_], Task, Res):- getPC(Schedule, Task, Res).

getPC(schedule(C, Ts), Task, process_cost(Task, C, ET)):-
  member(Task, Ts),
  process_cost(Task, C, ET), !.

getPC([_|Schedules], Task, Res):- getPC(Schedules, Task, Res), !.

getPC([], _, _):- fail.


% construct edges for a graph for a given schedule
make_io_edges(ScheduleList, Schedule, Res):-
  make_io_edges(ScheduleList, Schedule, [], Res).

make_io_edges(_, schedule(_, []), Edges, EdgesSet):-
  list_to_set(Edges, EdgesSet).

make_io_edges(Schedules, schedule(C, [T]), Edges0, Edges):-
  process_cost(T, C, ET),
  findall(process_cost(T, C, ET)-PC,
    make_edge(Schedules, [process_cost(T, C, ET)-PC]),
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
    make_edge(Schedules, [process_cost(T1, C, ET1)-PC]),
  OtherEdges),

  append( OtherEdges, Edges1, Edges2),
  make_io_edges(Schedules, schedule(C, [T2|Ts]), Edges2, Res), !.


% create graph for a schedule, uses SWI prolog's digraph implementation
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


% get distance or temporary costs from a tuple list
get_dist(PC, [tuple(PC, W)], W):- !.
get_dist(PC, [tuple(PC, W)|_], W):- !.

get_dist(PC, [tuple(_, _)|Rest], R):-
  get_dist(PC, Rest, R).

get_dist(_, [], _):- !.

% calculationg communication costs
get_comm_cost(process_cost(CurT, CurC, _), process_cost(NbT, NbC, _), Comm):-
  (CurC\=NbC, (channel(CurC, NbC, Lat, Bwidth), depends_on(NbT, CurT, Data)) ->
      Comm = Lat + Data/Bwidth
    ;
      Comm = 0
  ).


% add distance to a node to the distance list
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


% get longest path in graph realized with the distance list measured by ET
% writes the ETs for tasks TS in NewLenTo
% here LenTo can be reused form previous run in MinETList in find_heuristically
% saves lots of time
my_longest_path(G, Ts, LenTo, NewLenTo, Max):-
  longest_path(G, Ts, LenTo, NewLenTo),
  maplist(extractET, NewLenTo, ETList),
  max_list(ETList, Max).

% standard version of longest_path, for more info see:
% https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs_and_critical_paths
longest_path(G, Max):-
  top_sort(G, PCs),
  maplist(init_length_to, PCs, LengthTo),
  longest_path(G, PCs, LengthTo, List),
  maplist(extractET, List, ETList),
  max_list(ETList, Max).

longest_path(_, [], LengthTo, LengthTo):- !.

longest_path(G, [Cur|Rest], LenTo, R) :-
  get_dist(Cur, LenTo, CurW),
  children(G, Cur, Nbs),
  add_dist(tuple(Cur, CurW), Nbs, LenTo, LenTo, [], NewLenTo),
  longest_path(G, Rest, NewLenTo, R).


% predicate for tests
test(ET):-
  dep_sort_tasks(SortedTasks),
  add_pc_list_to_schedule(SortedTasks, [], [], ScheduleList),
  execution_time(ScheduleList, ET).
  % find_heuristically(ScheduleList),
