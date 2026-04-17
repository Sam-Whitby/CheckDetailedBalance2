(* ================================================================
   2D Reptation VMMC on a periodic square lattice
   ================================================================

   State:   flat array of length L^2 (row-major).
            State[[s]] = 0 (empty) or k ∈ {1,...,N} (labeled particle).
            L is inferred as Sqrt[Length[state]]; only IDs whose decoded
            array has perfect-square length are accepted by BitsToState.

   Site numbering (row-major, 1-indexed):
     1  2  ... L
     L+1 ...  2L
     ...
     (L-1)L+1 ... L^2

   Move:    Choose a random occupied site as the seed particle. Propose
            a virtual displacement of one lattice unit in a randomly
            chosen direction d ∈ {right, left, down, up}.

            If the target site is occupied, the move is rejected
            immediately (no hard-sphere seed moves).

            Otherwise the seed moves to the empty target, creating a
            vacancy at its old position. A reptation chain is then grown
            depth-first: at each vacancy the set of adjacent non-cluster
            occupied particles is found; one is chosen uniformly at
            random and moves into the vacancy (becoming a cluster
            member), creating a new vacancy at its own old position.
            Growth continues until the current vacancy has no adjacent
            non-cluster particles.

            After the chain is built a standard Metropolis acceptance
            step is applied on the total energy change of the proposed
            configuration.

   Reptation: each particle in the chain moves exactly one lattice unit
            into the space vacated by the previous particle, producing
            snake-like collective motion. A chain of length one is a
            plain single-particle translation; a chain of length two is
            a nearest-neighbour Kawasaki-style swap between two vacated
            and filled sites.

   Nearest-neighbour coupling:
            Same J_ab symbols as kawasaki_2d.wl / vmmc_2d.wl, summed
            over all horizontal and vertical bonds.

   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_2d.wl) ----------- *)

$cL[L_]        := $cL[L]     = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]  = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_] := $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]


(* ---- 2D lattice helpers -------------------------------------------------- *)

$right2D[s_, L_] := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c,L] + 1]
$down2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r,L]*L + c]
$left2D[s_, L_]  := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, (r-1)*L + Mod[c-2,L] + 1]
$up2D[s_, L_]    := With[{r = Ceiling[s/L], c = Mod[s-1,L]+1}, Mod[r-2,L]*L + c]

$allNeighbors2D[s_, L_] := {$right2D[s,L], $left2D[s,L], $down2D[s,L], $up2D[s,L]}

$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1


(* ---- Coupling function --------------------------------------------------- *)

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy -------------------------------------------------------------- *)

$uniqueBonds2D[L_] := $uniqueBonds2D[L] =
  DeleteDuplicates[Sort /@ Flatten[
    Table[{{s, $right2D[s,L]}, {s, $down2D[s,L]}}, {s, L^2}], 1]]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[couplingJ[state[[#[[1]]]], state[[#[[2]]]], 1] & /@ $uniqueBonds2D[L]]]


(* ---- Reptation chain builder --------------------------------------------- *)

(* Grows a reptation chain via DFS sequential vacancy-filling.
   seed moves to the empty adjacent site target = seed + dir.
   Returns {cluster, moveMap} where moveMap is an Association
   from each cluster site to its proposed new site.
   Returns None if the target site is occupied (hard-sphere rejection). *)

$reptationBuildChain[state_, L_, seed_, dir_] :=
  Module[{
    target    = $applyDir[seed, dir, L],
    cluster, inCluster, moveMap,
    vacancy, candidates, chosen
  },
    If[state[[target]] =!= 0, Return[None]];

    cluster   = {seed};
    inCluster = <|seed -> True|>;
    moveMap   = <|seed -> target|>;
    vacancy   = seed;

    While[True,
      (* Non-cluster occupied neighbours of the current vacancy *)
      candidates = DeleteDuplicates @ Select[
        $allNeighbors2D[vacancy, L],
        state[[#]] =!= 0 && !KeyExistsQ[inCluster, #] &
      ];
      If[candidates === {}, Break[]];

      (* Choose one uniformly at random to fill the vacancy *)
      chosen            = RandomChoice[candidates];
      AppendTo[cluster,   chosen];
      inCluster[chosen] = True;
      moveMap[chosen]   = vacancy;   (* chosen moves into the vacancy *)
      vacancy           = chosen     (* new vacancy = chosen's old site *)
    ];

    {cluster, moveMap}
  ]


(* ---- Algorithm ----------------------------------------------------------- *)

Algorithm[state_List] :=
  Module[{
    L, occupied, seed, dirs, dir,
    result, cluster, moveMap, newState, dE
  },
    L        = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];
    dirs = {{0,1},{0,-1},{1,0},{-1,0}};
    dir  = RandomChoice[dirs];

    result = $reptationBuildChain[state, L, seed, dir];
    If[result === None, Return[state]];   (* target occupied: hard-sphere reject *)

    {cluster, moveMap} = result;

    (* Apply moves: clear old sites first, then write new sites *)
    newState = state;
    Do[newState[[ cluster[[i]] ]] = 0,                                  {i, Length[cluster]}];
    Do[newState[[ moveMap[cluster[[i]]] ]] = state[[ cluster[[i]] ]],   {i, Length[cluster]}];

    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]
  ]


(* ---- Checker interface (identical to vmmc_2d.wl) ------------------------ *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]

numBeta = 1
