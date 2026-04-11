(* ================================================================
   2D rigid cluster move on a torus -- expected to FAIL
   ================================================================

   State:   flat L×L array (row-major, 1-indexed).
            State[[s]] = 0 (empty) or k (labeled particle).
            L inferred as Sqrt[Length[state]].

   Move:    Find all 4-connected clusters of particles.
            Pick one cluster uniformly at random.
            Pick a direction (E/S/W/N) uniformly at random.
            If all sites the cluster would shift into are currently
            empty, move the whole cluster one step in that direction
            and apply Metropolis acceptance; otherwise stay put.

   Why it fails: when a cluster moves adjacent to another cluster the
   two merge, reducing the number of available clusters in the new
   state. The proposal probability (1/n_clusters) then differs between
   the forward transition and its reverse, violating detailed balance.

   Concretely: state i has n_c(i) clusters, state j has n_c(j). If
   n_c(i) != n_c(j) then T(i→j) * π(i) != T(j→i) * π(j) in general.
   The symbolic checker detects this for the small systems tested by
   MaxBitString=1111111 (L=2, 2x2 torus).
   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_2d.wl) ----------- *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

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

(* Shift site s one step in direction dir (0=E, 1=S, 2=W, 3=N) *)
$shiftSite2D[s_, dir_, L_] :=
  Which[dir == 0, $right2D[s, L],
        dir == 1, $down2D[s, L],
        dir == 2, $left2D[s, L],
        True,     $up2D[s, L]]


(* ---- Coupling function --------------------------------------------------- *)

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[Table[
      couplingJ[state[[s]], state[[$right2D[s, L]]], 1] +
      couplingJ[state[[s]], state[[$down2D[s, L]]], 1],
      {s, Length[state]}]]]


(* ---- Cluster finder (4-connected BFS) ------------------------------------ *)

$findClusters2D[state_List, L_Integer] :=
  Module[{n = L*L, visited = ConstantArray[False, L*L],
          allClusters = {}, c, queue, s},
    Do[
      If[state[[seed]] != 0 && !visited[[seed]],
        c = {}; queue = {seed};
        While[queue != {},
          s = First[queue]; queue = Rest[queue];
          If[visited[[s]], Continue[]];
          visited[[s]] = True; AppendTo[c, s];
          Do[If[state[[nb]] != 0 && !visited[[nb]], AppendTo[queue, nb]],
             {nb, {$right2D[s,L], $down2D[s,L], $left2D[s,L], $up2D[s,L]}}]];
        AppendTo[allClusters, Sort[c]]],
      {seed, Range[n]}];
    allClusters]


(* ---- Algorithm ----------------------------------------------------------- *)

(* Pick a cluster uniformly (prob 1/n_clusters), pick a direction uniformly
   (prob 1/4), move if all target sites are empty, apply Metropolis.

   Fails DB: when a cluster moves adjacent to another cluster they merge,
   changing n_clusters between states i and j, so the 1/n_clusters proposal
   factors do not cancel in T(i→j)*π(i) = T(j→i)*π(j). *)
Algorithm[state_List] :=
  Module[{L, clusters, nC, cIdx, cluster, dir, targets, newFree, newState, dE},
    L        = Round[Sqrt[Length[state]]];
    clusters = $findClusters2D[state, L];
    nC       = Length[clusters];
    If[nC == 0, Return[state]];
    cIdx    = RandomInteger[{0, nC - 1}];
    cluster = clusters[[cIdx + 1]];
    dir     = RandomInteger[{0, 3}];
    targets = Map[$shiftSite2D[#, dir, L] &, cluster];
    (* Sites the cluster moves INTO that are not already in the cluster *)
    newFree = Complement[targets, cluster];
    If[AnyTrue[newFree, state[[#]] != 0 &], Return[state]];   (* blocked *)
    newState = state;
    Do[newState[[s]] = 0, {s, cluster}];   (* clear old positions *)
    Do[newState[[targets[[i]]]] = state[[cluster[[i]]]], {i, Length[cluster]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

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

numBeta = 1
