(* ================================================================
   1D rightward cluster move -- FAILS detailed balance
   ================================================================

   State:   flat array of length L (periodic ring).
            State[[i]] = 0 (empty) or k (labeled particle).

   Move:    Find all maximal connected runs of particles.
            Pick one uniformly at random.
            If the site immediately to its right is empty, slide
            the entire cluster one step clockwise and apply Metropolis.

   Why it fails: the reverse move (sliding left) is never attempted.
   For any transition i→j caused by a rightward slide,
   T(j→i) = 0 while T(i→j) > 0. The symbolic checker detects this
   as a violation for any component containing such a transition.
   ================================================================ *)


(* ---- Bijective integer encoding (identical to kawasaki_1d.wl) ----------- *)

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
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Cluster helpers ----------------------------------------------------- *)

(* Find all maximal connected runs of non-zero sites on an L-site ring.
   Returns a list of site lists, e.g. {{1,2},{5}} for a 6-site ring.
   Handles periodic boundary: sites L and 1 are adjacent, so if both
   are occupied the two runs are merged into one wrap-around cluster.
   Wrap-around clusters are stored as {last_run_sites..., first_run_sites...},
   so Last[cluster] is the clockwise endpoint and
   Mod[Last[cluster], L] + 1 is the next (rightmost) site. *)
$findClusters1D[state_List] :=
  Module[{L = Length[state], runs = {}, start = 0, inRun = False},
    Do[
      Which[
        state[[i]] != 0 && !inRun, start = i; inRun = True,
        state[[i]] == 0 && inRun,  AppendTo[runs, Range[start, i-1]]; inRun = False],
      {i, 1, L}];
    If[inRun, AppendTo[runs, Range[start, L]]];   (* close run ending at L *)
    (* Merge wrap-around: first run starts at site 1 AND last run ends at site L *)
    If[Length[runs] >= 2 && runs[[1,1]] == 1 && runs[[-1,-1]] == L,
      Join[{Join[runs[[-1]], runs[[1]]]}, runs[[2 ;; -2]]],
      runs]]

(* Slide cluster one step clockwise (each site s → Mod[s,L]+1).
   Saves particle values before clearing old sites to handle any wrap-around. *)
$slideRight1D[state_List, cluster_List] :=
  Module[{L = Length[state], newState = state, vals},
    vals = state[[cluster]];
    Do[newState[[s]] = 0, {s, cluster}];
    Do[newState[[Mod[cluster[[i]], L] + 1]] = vals[[i]], {i, Length[cluster]}];
    newState]


(* ---- Algorithm ----------------------------------------------------------- *)

(* Pick a cluster uniformly; if its right neighbour is empty, slide it right.
   Never proposes a leftward slide → breaks detailed balance. *)
Algorithm[state_List] :=
  Module[{L, clusters, nC, b, cluster, rightSite, newState, dE},
    L        = Length[state];
    clusters = $findClusters1D[state];
    nC       = Length[clusters];
    If[nC == 0, Return[state]];
    b         = RandomInteger[{0, nC - 1}];
    cluster   = clusters[[b + 1]];
    rightSite = Mod[Last[cluster], L] + 1;   (* site clockwise of cluster *)
    If[state[[rightSite]] != 0, Return[state]];   (* blocked *)
    newState = $slideRight1D[state, cluster];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
